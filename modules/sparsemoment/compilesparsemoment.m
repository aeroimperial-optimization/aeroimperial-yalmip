function prog = compilesparsemoment(x, p, h, g, omega, mass, options, cliques)

% Compile conic program corresponding to a sparsity-exploiting moment-SOS
% relaxation of a POP. The conic problem is in the form
%
% min   -b'*y
% s.t.  c - At*y \in K
%
% where K is a cone. For moment-SOS relaxations, the cones are:
%    * K.f: the zero cone (of dimension equal to K.f)
%    * K.l: the nonnegative orthant (of dimension equal to K.l)
%    * K.s: a collection of semidefinite cones (each of linear size K.s(i))

% ============================================================================ %
%                               INITIALIZATION                                 %
% ============================================================================ %
% Inputs correct?
if nargin < 5; omega = []; end                  % No relaxation order? Set as empty
if nargin < 6; mass = 1; end                    % No mass specified? Use 1
if nargin < 7; options = sdpsettings; end       % No options? use yalmip default
if nargin < 8; cliques = []; end                % No cliques? set to empty
if isempty(options); options = sdpsettings; end % Empty options? use yalmip default
if isempty(mass); mass = 1; end                 % Empty mass? set to 1

% Display a nice header if needed
if options.verbose
    disp(repmat('=',1,40))
    disp('    MOMENT RELAXATION WITH SPARSITY')
    disp(repmat('=',1,40))
end

% Start the clock
starttime = tic;

% Get number of variables in POP and yalmip internal IDs
numx = numel(x);
xID = depends(x);

% Check if omega is large enough, otherwise increase it
d = max(degree([h;g]), degree(p));
d = ceil(d/2);
if isempty(omega) || (omega < d)
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
    CD = corrSparsityCliques(x, p, [h; g]);
else
    CD.Set = cliques;
    CD.NoC = numel(CD.Set);
    clear cliques
end


% List all constraints
cnstr = [h; g];
num_h = numel(h);
num_g = numel(g);
num_cnstr = num_h + num_g;

% Create an empty model for the problem
% The vector all_moments lists all moments in the relaxation
% The constraints are in sedumi format for simplicity: c - At*y \in K
[all_moments,At,c,K,momentID,gramMonomials] = initializeModel(numx,omega,CD);
isMomentMatrix = cell(CD.NoC,1);
hash = rand(numx,1);
all_moments_hash = all_moments*hash;
[all_moments_hash, idx] = uniquetol(all_moments_hash, 1e-12);
all_moments = all_moments(idx, :);
[all_moments, sortidx] = sortrows(all_moments); % ensure consistency between runs
all_moments_hash = all_moments_hash(sortidx);
num_moments = size(all_moments, 1);

% ============================================================================ %
%                            BUILD OBJECTIVE                                   %
% ============================================================================ %
% Build objective as -b'*y - b0, where y is the vector of moments except for
% the zero moment. The constant term is ignored and added at the end
if options.verbose
    disp('Constructing the objective function...')
end
[powers, b] = getexponentbase(p,x);
[ia,ib] = ismembertol(powers*hash, all_moments_hash);
b0 = -full(b(~ia));
if isempty(b0); b0 = 0; end
b = sparse(ib(ia), 1, -b(ia), num_moments, 1);


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
for i = 1:CD.NoC
    
    % Display progress if required
    if options.verbose
        if i==1
            fprintf('Constructing the constraints (%i cliques)...\n',CD.NoC)
            fprintf('Progress: ')
        else
            fprintf(repmat('\b',1,10))
        end
        
        fprintf(repmat('%%',1,floor(10*i/CD.NoC)))
        fprintf(repmat(' ',1,10-floor(10*i/CD.NoC)))
    end

    % Variable ID and symbolic variables in this clique
    % Use low-level YALMIP command for fast construction
    var_id = xID(CD.Set{i});
    num_vars = length(var_id);
    var_base = spdiags(ones(num_vars,1),1,num_vars,num_vars+1);
    xloc = sdpvar(num_vars, 1, [], var_id, var_base);
    
    % Find moments that belong to this clique
    idx = setdiff(1:numx,CD.Set{i});
    inclique{i} = find(sum(all_moments(:,idx),2)==0);
    clique_moments = all_moments(inclique{i},:);
    
    % Get the basic moment LMI without using symbolic variables
    % Code is hard to read, but essentially work with powers of monomials
    % and reduce to random vector for speedy searches
    [At(i),c(i),K(i),momentID{i},gramMonomials{i}] = buildMomentMatrix(num_vars,CD.Set{i},omega,numx,num_moments,clique_moments,inclique{i},[],At(i),c(i),K(i),mass);
    isMomentMatrix{i}(end+1) = 1;
    
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
            [At(i),c(i),K(i)] = buildConstraintMatrix(num_vars,CD.Set{i},omega,numx,num_moments,clique_moments,inclique{i},cdata,At(i),c(i),K(i),mass);
        else
            % Moment localizing matrices
            [At(i),c(i),K(i)] = buildMomentMatrix(num_vars,CD.Set{i},omega,numx,num_moments,clique_moments,inclique{i},cdata,At(i),c(i),K(i),mass);
            isMomentMatrix{i}(end+1) = 0;
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
    prog = makeConicUniqueMoments(At, b, c, K, isMomentMatrix, options);
    prog.CD = CD;
    prog.b0 = b0;
    prog.all_moments = all_moments;
    prog.gramMonomials = gramMonomials;
else
    % TO DO: use matching variable
    if options.verbose
        fprintf('\nSetting up conic problem (use splitting variable)...'); 
    end
    prog = makeConicSplittingVariable(At, b, c, K, isMomentMatrix, inclique, num_moments, options);
    prog.CD = CD;
    prog.b0 = b0;
    prog.all_moments = all_moments;
    prog.gramMonomials = gramMonomials;
end

% Clock setup time
prog.setuptime = toc(starttime);
fprintf('\nCompilation complete!\n'); 