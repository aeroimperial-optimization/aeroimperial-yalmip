function [pstar, y, unique_exponents, sol, model] = solvesparsemoment(x, p, h, g, omega, MassValue, options, cliques)

% [pstar,y,exponents,sol,model] = solvesparsemoment(x,p,h,g,omega,maxMass,options,cliques)
%       solves the moment relaxation of the standard-form polynomial 
%       optimization problem (POP)
%
%       min_x p(x) s.t. h_i(x)=0, g_j(x)>=0
%
%       by exploiting correlative sparsity. 
%
%       NOTE: This is research code. It is not optimized, not guaranteed to
%       work, and probably buggy. Use at your own risk!
%
%       Required inputs
%       ----------------
%       * x: the optimization variables of the POP (sdpvar vector)
%       * p: the polynomial objective to be MINIMIZED
%       * h: a vector of polynomial equality constraints, h_i(x)=0
%       * g: a vector of polynomial inequality constraints, g_i(x)>=0
%       * omega: the degree of the relaxation. The max degree of monomials
%                considered in the moment relaxation and associated SOS
%                program will be 2*omega.
%
%       Optional inputs
%       ---------------
%       * MassValue: total mass of the the moment-generating measure. The default
%                    value is 1, but setting a different positive value may
%                    improve the numerical behaviour
%       * options: yalmip options created with sdpsettings()
%       * cliques: sets of variables specifying the problem's sparsity
%                  structure. If not provided, the sparsity is detected
%                  automatically using a chordal extension of the
%                  correlative sparsity graph (see Waki et al, 2006)
%
%       Outputs
%       -------
%       * pstar: optimal value of the moment-SDP relaxation. it is a lower
%                bound on the optimal value of the POP
%       * y: the vector of moments solving the moment-SDP relaxation
%       * sol: the solver's output (if required)
%       * model: the conic problem model (in SeDuMi format)
%       * exponents: the exponents of the monomials modelled by y in the
%                    moment-SDP relaxation. For POPs with a unique solution
%                    and whose sparsity satisfies the running intersection
%                    property, one has
%               
%                    y(i) -> x(1)^exponents(i,1) * ... * x(n)^exponents(i,n)
%
%                    as the relaxation parameter omega tends to infinity.
%       
%
% Giovanni Fantuzzi
% 02 Dec 2019

% Optional inputs
if nargin < 6; MassValue = 1; end             % No max mass bound? use 1
if nargin < 7; options = sdpsettings; end   % No options? use yalmip default
if nargin < 8; cliques = []; end            % No cliques? set to empty
if isempty(options); options = sdpsettings; end % Empty options? use yalmip default


if options.verbose
    disp(repmat('=',1,40))
    disp('    MOMENT RELAXATION WITH SPARSITY')
    disp(repmat('=',1,40))
end

% Clock
starttime = tic;

% get number of variables in POP
numx = numel(x);

% Check if omega is large enough, otherwise increase it
d = max(degree([h;g]), degree(p));
d = ceil(d/2);
if omega < d
    omega = d;
    if options.verbose
        fprintf('Specified relaxation parameter omega too small. Using omega = %i\n',omega);
    end
end

% Get cliques of the correlative sparsity pattern if not provided
% The current code uses a chordal extension of the correlative sparsity
% pattern, which may give very large cliques!
if isempty(cliques)
    if options.verbose
        disp('Detecting correlative sparsity...');
    end
    cliques = corrSparsityCliques(x, p, [h; g]);
end
num_cliques = numel(cliques);

% list all constraints
cnstr = [h; g];
num_h = numel(h);
num_g = numel(g);
num_cnstr = num_h + num_g;

% Get moment indices (= exponents of monomials) and base of the objective
% Objective will be written as b'*y eventually, here get a temporary b to
% be inflated to include all moments in the constraints later
if options.verbose
    disp('Constructing the objective function...')
end
[exponents_y, b] = getexponentbase(p,x);
b = b(:);

% Initialize cone K (sedumi format)
% Include various cones that YALMIP can handle even if unused
K.f = 0;                    % free variables (to be populated)
K.l = 0;                    % nonnegative variables (to be populated)
K.s = [];                   % semidefinite variables (to be populated)
K.q = 0;                    % to be ignored
K.e = 0;                    % to be ignored
K.c = 0;                    % to be ignored
K.r = 0;                    % to be ignored
K.p = 0;                    % to be ignored
K.m = 0;                    % to be ignored
K.scomplex = [];            % to be ignored
K.xcomplex = [];            % to be ignored 
K.sos = [];                 % to be ignored
K.schur_funs = [];          % to be ignored
K.schur_data = [];          % to be ignored
K.schur_variables = [];     % to be ignored

% Handle constraints: loop over cliques to exploit sparsity, and loop over
% constraints to decide which ones correspond to a given clique
if options.verbose
    fprintf('Constructing the constraints (%i cliques)...\n',num_cliques)
end
At.f = {};       % constraints c - At*y==0
At.l = {};       % constraints c - At*y \in K.l
At.s = {};       % constraints c - At*y \in K.s
c.f = [];
c.l = [];
c.s = [];
exponents_clique = [];
exponents_clique_lmi = [];
whichClique.f = [];
whichClique.s = [];

% Find variables in each constrains and degree of constraints
constr_vars = cell(num_cnstr,1);
constr_degs = zeros(num_cnstr,1);
for i = 1:num_cnstr
    constr_vars{i} = depends(cnstr(i));
    constr_degs(i) = degree(cnstr(i));
end

% Loop over cliques
for i = 1:num_cliques
    
    % Display progress?
    if options.verbose
        if i==1
            fprintf('Progress: ')
        else
            fprintf(repmat('\b',1,10))
        end
        
        fprintf(repmat('%%',1,floor(10*i/num_cliques)))
        fprintf(repmat(' ',1,10-floor(10*i/num_cliques)))
    end

    % ID and symbolic variables in this clique
    % Also construct monomials in this clique up to degree omega and
    % PSD matrix of monomials (needed to build moment localizing matrices)
    var_id = cliques{i};
    num_vars = length(var_id);
    var_base = spdiags(ones(num_vars,1),1,num_vars,num_vars+1);
    xloc = sdpvar(num_vars, 1, [], var_id, var_base);
    
    % -------------------------------------------------------
    % Get the basic moment LMI without using symbolic variables
    % Code is hard to read, but essentially work with powers of monomials
    % and reduce to random vector for speedy searches
    Mpow = monolistcoeff(num_vars,omega,omega);         % get exponents for local vars
    K.s(end+1) = size(Mpow,1);                          % size of cone (linear)
    [Mpow2,N_unique] = monomialproducts({Mpow});
    [ii,jj,vv] = find(N_unique(:,3:end));                   % find
    clique_all_exponents{i} = sparse(ii,var_id(jj),vv,max(ii),numx); % inflate to exponents of all vars
    hash = randn(num_vars,1);
    Nhash = N_unique(:,3:end)*hash;     % vectorize with randomness for speed
    Mpow2hash = Mpow2{1}(:,3:end)*hash;   % faster to search for matching numbers
    rows = [];
    cols = [];
    Mpow2 = Mpow2{1};
    [~,TEMP] = ismembertol(Mpow2hash, Nhash);
    rows = [rows; Mpow2(:,1)+(Mpow2(:,2)-1).*K.s(end)];
    cols = [cols; TEMP];
    c.s = [c.s; sparse(K.s(end)^2, 1)];
    At.s{end+1} = sparse(rows,cols,-1,K.s(end)^2,length(Nhash));
    whichClique.s = [whichClique.s, i];
    % -------------------------------------------------------
    
    % Find which constraints belong to this clique
    cnstr_indices = cellfun(@(C)fastismember(C, var_id), constr_vars);
    cnstr_in_clique = find(cnstr_indices);
    num_cnstr_in_clique = length(cnstr_in_clique);
    
    % Add indices of cliques to list
    % whichClique = [whichClique; ones(num_cnstr_in_clique,1).*i];
    
    % Loop over constraints and build equalities or moment localizing LMIs
    % if constraint h_j or g_j is linked to the current clique. Constraints
    % are represented as c-At*y \in K* (the dual of K). Here, y has
    % repeated entries that will be eliminated later!
    for kindex = 1:num_cnstr_in_clique
        % Get value of j
        j = cnstr_in_clique(kindex);
        % Bingo! Proceed depending on whether it is an equality or
        % inequality constraint
        if j <= num_h
            % HARD TO READ BUT AVOID EXPENSIVE SYMBOLIC OPERATIONS
            % Essentially, multiply the list of monomials by the constraint
            % and derive equalities. Work with powers and use random hasing
            % for speedy searches
            degz = 2*omega-constr_degs(j);
            zpow = monolistcoeff(num_vars,degz,degz);
            [h_pow, h_coef] = getexponentbase(cnstr(j),xloc);
            rows = [];
            cols = [];
            vals = [];
            for kk = 1:length(h_coef)
                Q = zpow + repmat(h_pow(kk,:),size(zpow,1),1);
                Qhash = Q*hash;
                rows = [rows; (1:length(Qhash))'];
                vals = [vals; -h_coef(kk).*ones(length(Qhash),1)];
                [~,TEMP] = ismembertol(Qhash, Nhash);
                cols = [cols; TEMP];
            end
            c.f = [c.f; sparse(length(Qhash), 1)];
            At.f{end+1} = sparse(rows,cols,vals,length(Qhash),length(Nhash));
            whichClique.f = [whichClique.f, i];
            K.f = K.f + length(Qhash);
            % ---------------------------------------------------
        else
            % AVOID EXPENSIVE SYMBOLIC OPERATIONS
            % Essentially, multiply the matrix of monomials Q by the constraint
            % and derive the moment localizing matrix. 
            % Work with powers and use random hasing for speedy searches
            degQ = omega - ceil( 0.5*constr_degs(j) );
            Qpow = monolistcoeff(num_vars,degQ,degQ);
            nsdp = size(Qpow,1);
            [Qpow,Q_unique] = monomialproducts({Qpow});
            [g_pow, g_coef] = getexponentbase(cnstr(j),xloc);
            K.s(end+1) = nsdp;
            rows = [];
            cols = [];
            vals = [];
            new_exponents = [];
            for kk = 1:length(g_coef)
                Q = Qpow{1} + repmat([0, 0, g_pow(kk,:)],size(Qpow{1},1),1);
                P = Q_unique + repmat([0, 0, g_pow(kk,:)],size(Q_unique,1),1);
                Phash = P(:,3:end)*hash;   % vectorize with randomness for speed
                Qhash = Q(:,3:end)*hash;   % faster to search for matching numbers
                [~,LOCB1] = ismembertol(Qhash, Phash);
                [~,LOCB2] = ismembertol(Phash, Nhash);
                TEMP = LOCB2(LOCB1);
                rows = [rows; Q(:,1)+(Q(:,2)-1).*nsdp];
                cols = [cols; TEMP];
                vals = [vals; -g_coef(kk).*ones(size(TEMP,1),1)];
            end
            c.s = [c.s; sparse(nsdp^2, 1)];
            At.s{end+1} = sparse(rows,cols,vals,nsdp^2,length(Nhash));
            whichClique.s = [whichClique.s, i];
            % ---------------------------------------------------
        end % End if equality/inequality
    end % end for loop over all constraints
end % end loop over cliques

% Add constraint on mass
% MassValue - y0 == 0
K.f = K.f + 1;
c.f = [c.f; MassValue];
At.f{end+1} = 1;

% Now combine moments from constraints and objective and strip duplicate
% monomials to construct conic problem. Overwrite exponents_y since not
% needed from now on.
shift = size(exponents_y,1);
exponents_y = [{exponents_y}, {sparse(1,length(x))}, clique_all_exponents];
if options.sparsemoment.mergeCliques
    % Remove duplicate copies of the moments and build the constraints
    % Transform At from cell array into single sparse matrix
    if options.verbose
        fprintf('\nSetting up conic problem (removing duplicate moments)...'); 
    end
    [At, b, c, K, unique_exponents] = makeConicUniqueMomentsNew(At, b, c, K, exponents_y, shift, whichClique);
else
    % TO DO: use matching variable
    if options.verbose
        fprintf('\nSetting up conic problem (use splitting variable)...'); 
    end
    [At, b, c, K, unique_exponents] = makeConicSplittingVariable(At, b, c, K, exponents_y, shift, whichClique, clique_all_exponents);
end

% Clock setup time
setuptime = toc(starttime);

% Finally, solve and set outputs. Try to use generic yalmip format to 
% interface easily with user-preferred choice. Tested with:
% sedumi, mosek, sdpt3, cdcs, scs 
% Could be buggy with other solvers.
% NOTE: ensure we save solver outputs to get the solution
if options.verbose; disp('Solving!'); end
options.savesolveroutput = 1;
[solver,problemClass] = sparsemoments_getsolvers(options);
interfacedata = sparsemoments_interfacedata(At,b,c,K,options,solver,problemClass);
try
    eval(['output = ' solver.call '(interfacedata);']);
catch
    error('Woops, something went wrong in the solver!')
end

% Get solution and scale y by the zero-th moment
y = output.Primal;
y = y./MassValue;

% Get optimal value of the moments-SDP relaxation
% Need to use y with right scaling here. We minimize b'*y (difference
% between minimization and maximization is taken into account previously by
% handling minus signs)
pstar = b.'*y;

% Output solver solution if desired
if nargout > 2
    sol = output.solveroutput; 
    sol.setupTime = setuptime;
end

% Output model if needed
% Change sign of b here to have dual-standard-form problem:
% max b'*y    s.t.    c - At*y \in K
if nargout > 3
    model.At = At;
    model.b = -b;
    model.c = c;
    model.K = K;
end