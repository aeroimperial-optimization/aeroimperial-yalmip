function [pstar, y, all_moments, sol, model] = solvesparsemoment(x, p, h, g, omega, mass, options, cliques)

% [pstar,y,exponents,sol,model] = solvesparsemoment(x,p,h,g,omega,mass,options,cliques)
%       solves the moment relaxation of the standard-form polynomial 
%       optimization problem (POP)
%
%       min_x p(x) s.t. h_i(x)=0, g_j(x)>=0
%
%       by exploiting correlative sparsity. 
%
%       NOTES: 
%       (1) This is research code. It is not optimized, not guaranteed to
%           work, and probably buggy. Use at your own risk!
%       (2) This code does NOT exploit symmetry
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
%                    improve the numerical behaviour (it amounts to problem
%                    scaling)
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

% ============================================================================ %
%                               INITIALIZATION                                 %
% ============================================================================ %
% Optional inputs
if nargin < 6; mass = 1; end                    % No mass specified? Use 1
if nargin < 7; options = sdpsettings; end       % No options? use yalmip default
if nargin < 8; cliques = []; end                % No cliques? set to empty
if isempty(options); options = sdpsettings; end % Empty options? use yalmip default

% Empty mass?
if isempty(mass); mass = 1; end

% Display a nice header if needed
if options.verbose
    disp(repmat('=',1,40))
    disp('    MOMENT RELAXATION WITH SPARSITY')
    disp(repmat('=',1,40))
end

% Start the clock
starttime = tic;

% Get number of variables in POP
numx = numel(x);

% Check if omega is large enough, otherwise increase it
d = max(degree([h;g]), degree(p));
d = ceil(d/2);
if omega < d
    omega = d;
    if options.verbose
        fprintf('Specified relaxation parameter omega too small (POP has degree %i)!\nUsing omega = %i\n',2*d,omega);
    end
end

% Get cliques of the correlative sparsity pattern if not provided
% NOTE: the clique detection uses a standard chordal extension of the 
% correlative sparsity graph, which may give very large cliques.
if isempty(cliques)
    if options.verbose
        disp('Detecting correlative sparsity...');
    end
    CD.cliques = corrSparsityCliques(x, p, [h; g]);
else
    CD.cliques = cliques;
    clear cliques
end
CD.num_cliques = numel(CD.cliques);

% List all constraints
cnstr = [h; g];
num_h = numel(h);
num_g = numel(g);
num_cnstr = num_h + num_g;

% Create an empty model for the problem
% The vector all_moments lists all moments in the relaxation
% The constraints are in sedumi format for simplicity: c - At*y \in K
[all_moments,At,c,K,momentID] = initializeModel(numx,omega,CD);
hash = rand(numx,1);
all_moments_hash = all_moments*hash;
[all_moments_hash, idx] = uniquetol(all_moments_hash, 1e-12);
all_moments = all_moments(idx, :);
num_moments = size(all_moments, 1);

% ============================================================================ %
%                            BUILD OBJECTIVE                                   %
% ============================================================================ %
% Build objective as b'*y + b0, where y is the vector of moments except for
% the zero moment. The constant term is ignored and added at the end
if options.verbose
    disp('Constructing the objective function...')
end
[powers, b] = getexponentbase(p,x);
[ia,ib] = ismembertol(powers*hash, all_moments_hash);
b0 = full(b(~ia));
if isempty(b0); b0 = 0; end
b = sparse(ib(ia), 1, b(ia), num_moments, 1);


% ============================================================================ %
%                            BUILD CONSTRAINTS                                 %
% ============================================================================ %
% Find variables in each constrains and the degree of constraints
constr_vars = cell(num_cnstr,1);
constr_degs = zeros(num_cnstr,1);
for i = 1:num_cnstr
    constr_vars{i} = depends(cnstr(i));
    constr_degs(i) = degree(cnstr(i));
end

% Actually build the constraints by looping over the cliques
for i = 1:CD.num_cliques
    
    % Display progress if required
    if options.verbose
        if i==1
            fprintf('Constructing the constraints (%i cliques)...\n',CD.num_cliques)
            fprintf('Progress: ')
        else
            fprintf(repmat('\b',1,10))
        end
        
        fprintf(repmat('%%',1,floor(10*i/CD.num_cliques)))
        fprintf(repmat(' ',1,10-floor(10*i/CD.num_cliques)))
    end

    % Variable ID and symbolic variables in this clique
    % Use low-level YALMIP command for fast construction
    var_id = CD.cliques{i};
    num_vars = length(var_id);
    var_base = spdiags(ones(num_vars,1),1,num_vars,num_vars+1);
    xloc = sdpvar(num_vars, 1, [], var_id, var_base);
    
    % Get the basic moment LMI without using symbolic variables
    % Code is hard to read, but essentially work with powers of monomials
    % and reduce to random vector for speedy searches
    [At(i),c(i),K(i),momentID{i}] = buildMomentMatrix(num_vars,var_id,omega,numx,all_moments,[],At(i),c(i),K(i),mass);
    
    % Find which constraints belong to this clique
    cnstr_indices = cellfun(@(C)fastismember(C, var_id), constr_vars);
    cnstr_in_clique = find(cnstr_indices);
    num_cnstr_in_clique = length(cnstr_in_clique);
    
    % Loop over constraints and build equalities or moment localizing LMIs
    % if constraint h_j or g_j is linked to the current clique. Constraints
    % are represented as c-At*y \in K.
    % NOTE: y has repeated entries that will be eliminated later
    for kindex = 1:num_cnstr_in_clique
        % Get value of j and constraint polynomial structure
        j = cnstr_in_clique(kindex);
        cdata.degree = constr_degs(j);
        [cdata.pows, cdata.coef] = getexponentbase(cnstr(j),xloc);
        if j <= num_h
            % Equality constraints
            [At(i),c(i),K(i)] = buildConstraintMatrix(num_vars,var_id,omega,numx,all_moments,cdata,At(i),c(i),K(i),mass);
        else
            % Moment localizing matrices
            [At(i),c(i),K(i)] = buildMomentMatrix(num_vars,var_id,omega,numx,all_moments,cdata,At(i),c(i),K(i),mass);
        end
        % End if equality/inequality
    end
    % end for loop over all constraints
end
% end loop over cliques


% ============================================================================ %
%                             COMPILE THE SDP                                  %
% ============================================================================ %
% Now combine moments from constraints and objective and strip duplicate
% monomials to construct conic problem. Overwrite exponents_y since not
% needed from now on.
if options.sparsemoment.mergeCliques
    % Remove duplicate copies of the moments and build the constraints
    % Transform At from cell array into single sparse matrix
    if options.verbose
        fprintf('\nSetting up conic problem (removing duplicate moments)...'); 
    end
    [At, b, c, K, PROJ, bshift, y0] = makeConicUniqueMoments(At, b, c, K, options);
else
    % TO DO: use matching variable
    error('No other option is known: use "mergeCliques" for now!')
    if options.verbose
        fprintf('\nSetting up conic problem (use splitting variable)...'); 
    end
    [At, b, c, K, PROJ, bshift, y0] = makeConicSplittingVariable(At, b, c, K, options);
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

% Get solution, scale y by the zero-th moment and get optimal value of the 
% moments-SDP relaxation.
y = y0 + PROJ*output.Primal;
y = y./mass;
pstar = (b.'*output.Primal + bshift)/mass + b0;

% Output solver solution if desired
if nargout > 3
    sol = output.solveroutput; 
    sol.setupTime = setuptime;
end

% Output model if needed
% Change sign of b here to have dual-standard-form problem:
% max b'*y    s.t.    c - At*y \in K
if nargout > 4
    model.At = At;
    model.b = -b;
    model.c = c;
    model.K = K;
end