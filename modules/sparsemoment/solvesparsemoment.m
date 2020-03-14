function [pstar, y, unique_exponents, sol, model] = solvesparsemoment(x, p, h, g, omega, maxMass, options, cliques)

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
%       * omega: the degree of the relaxation (equivalent to the degree of
%                the Lagrange multipliers for the S-procedure on the SOS side)
%
%       Optional inputs
%       ---------------
%       * maxMass: the maximum mass of the moment-generating
%                  measure. This makes the numerics behave a little 
%                  better compared to fixing the mass to be exactly equal
%                  to 1 (or any other constant). DEFAULT: 1
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
if nargin < 6; maxMass = 1; end             % No max mass bound? use 1
if nargin < 7; options = sdpsettings; end   % No options? use yalmip default
if nargin < 8; cliques = []; end            % No cliques? set to empty
if isempty(options); options = sdpsettings; end % Empty options? use yalmip default


if options.verbose
    disp(repmat('=',1,40))
    disp('    MOMENT RELAXATION WITH SPARSITY')
    disp(repmat('=',1,40))
end

% Odd omega?
if rem(omega,2) ~= 0
    omega = 2*ceil(0.5*omega);
    disp(['Input omega must be even. Increasing omega to ' num2str(omega)])
end

% Clock
starttime = tic;

% get number of variables in POP
numx = numel(x);

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

% Get total degree for general moment matrices; make even if odd
d = max(degree([h;g])+omega, degree(p));
d = ceil(d/2);

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
K.l = 2;                    % nonnegative variables (there will be two)
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
c.s = [];
exponents_clique = [];
exponents_clique_lmi = [];
whichClique.f = [];
whichClique.s = [];

% Find variables in each constrains
constr_vars = cell(num_cnstr,1);
for i = 1:num_cnstr
    constr_vars{i} = depends(cnstr(i));
end

% Loop over cliques
for i = 1:num_cliques
    
    %     % Display
    %     if options.verbose
    %         disp(['Clique number ',num2str(i)])
    %     end
    
    % ID and symbolic variables in this clique
    % Also construct monomials in this clique up to degree omega and
    % PSD matrix of monomials (needed to build moment localizing matrices)
    var_id = cliques{i};
    num_vars = length(var_id);
    var_base = spdiags(ones(num_vars,1),1,num_vars,num_vars+1);
    xloc = sdpvar(num_vars, 1, [], var_id, var_base);
    zpow = monolistcoeff(num_vars,omega,omega);
    Qpow = monolistcoeff(num_vars,omega/2,omega/2);
    nsdp = size(Qpow,1);
    [Qpow,Q_unique] = monomialproducts({Qpow});
    
    %     % Get the basic moment LMI
    %     % SYMBOLIC VERSION -- SLOW!
    %     % Remove the constant monomial, known (equivalent to solving y_0=1).
    %     M = monolist(xloc, d);
    %     M = M*M';
    %     [new_exponents, base] = getexponentbase(M, x);
    %     cs = [cs; sparse(size(base,1), 1)];
    %     Ats{end+1} = -base;
    %     exponents_clique_lmi = [exponents_clique_lmi; new_exponents];
    %     K.s(end+1) = size(M,1);
    
    % -------------------------------------------------------
    % Without using symbolic variables
    % Mpow = monpowers(num_vars,d);           % get exponents for local vars
    Mpow = monolistcoeff(num_vars,d,d);       % get exponents for local vars
    K.s(end+1) = size(Mpow,1);                % size of cone (linear)
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
    % exponents_clique_lmi = [exponents_clique_lmi; clique_all_exponents{i}];
    % -------------------------------------------------------
    
    % Find which constraints belong to this clique
%     cnstr_indices = cellfun(@(C)all(ismember(C, var_id)), constr_vars);
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
            %                 % Equality: get equality constraints on moments
            %                 % Do not remove redundant variables yet
            %                 temp = cnstr(j).*z; %%% SLOW -- CANNOT RESOLVE
            %                 [new_exponents, base] = getexponentbase(temp, x);
            %                 c = [c; sparse(size(base,1), 1)];
            %                 At{end+1} = -base;
            %                 exponents_clique = [exponents_clique; new_exponents];
            %                 K.f = K.f + size(base,1);
            % ---------------------------------------------------
            % AVOID EXPENSIVE SYMBOLIC OPERATIONS
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
            % exponents_clique = [exponents_clique; clique_all_exponents{i}];
            whichClique.f = [whichClique.f, i];
            K.f = K.f + length(Qhash);
            % ---------------------------------------------------
        else
            % Inequality: get PSD constraints on
            % Do not remove redundant variables yet
            %                 temp = cnstr(j).*Q; %%% SLOW -- CANNOT RESOLVE
            %                 [new_exponents, base] = getexponentbase(temp, x);
            %                 cs = [cs; sparse(size(base,1), 1)];
            %                 Ats{end+1} = -base;
            %                 exponents_clique_lmi = [exponents_clique_lmi; new_exponents];
            %                 K.s(end+1) = nsdp;
            % ---------------------------------------------------
            % AVOID EXPENSIVE SYMBOLIC OPERATIONS
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
            % exponents_clique_lmi = [exponents_clique_lmi; clique_all_exponents{i}];
            % ---------------------------------------------------
        end % End if equality/inequality
    end % end for loop over all constraints
end % end loop over cliques

% Add bounds on mass: two inequalities in cone K.l, modelled in standard
% conic form as: c-At*y \in K.l
% * mass is bounded: mass_bound - y_0 >= 0
% * mass in non-negative: y_0 >= 0
c.l = [maxMass; sparse(1,1)];
At.l{1} = 1;  % At for bound on mass
At.l{2} = -1; % At for non-negativity
%exponents_clique = [exponents_clique; sparse(2,length(x))];

% Now combine moments from constraints and objective and strip duplicate
% monomials to construct conic problem. Overwrite exponents_y since not
% needed from now on.
% -------------------------------------
% OLD VERSION
% exponents_clique = [exponents_clique; exponents_clique_lmi];
% [-b.'; sparse(size(exponents_clique,1),1)];
% At = [At, Ats];
% c = [c; cs];
% -------------------------------------
shift = size(exponents_y,1);
%exponents_y = [exponents_y; exponents_clique];
exponents_y = [{exponents_y}, {sparse(1,length(x))}, clique_all_exponents];
if options.sparsemoment.mergeCliques
    % Remove duplicate copies of the moments and build the constraints
    % Transform At from cell array into single sparse matrix
    if options.verbose
        disp('Setting up conic problem (removing duplicate moments)...'); 
    end
    [At, b, c, K, unique_exponents] = makeConicUniqueMomentsNew(At, b, c, K, exponents_y, shift, whichClique);
else
    % TO DO: use matching variable
    if options.verbose
        disp('Setting up conic problem (use splitting variable)...'); 
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
% ======================================================================= %
% OLD VERSION: CALL MOSEK OR CDCS DIRECTLY
% if isempty(options.solver); options.solver = 'mosek'; end
% switch options.solver
%     case 'mosek'
%         [y, sol] = mosekSDP(At,b,c,K,options);
%         sol.processing_time = setuptime;
%     case 'cdcs'
%         [~,y,sol] = cdcs(At,b,c,K,options.cdcs);
%     otherwise
%         error('Solver not recognized')
% end
% ======================================================================= %

% Get solution and scale y by the zero-th moment
y = output.Primal;
mass = y(sum(unique_exponents,2)==0);
y = y./mass;

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