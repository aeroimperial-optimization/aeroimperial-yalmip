function [all_moments,At,c,K,momentID] = initializeModel(numx,omega,cliqueData)
% Initalize an empty SeDuMi model for sparse moment relaxations

% Set stuff up
numRows = 0;
rows = [];
cols = [];
vals = [];
momentID = cell(cliqueData.num_cliques,1);
gramMonomials = cell(cliqueData.num_cliques,1);
% Loop over cliques
for i = 1:cliqueData.num_cliques
    
    % First, empty cones
    % Include all fields for compatibility, even though most will be ignored.
    K(i).f = 0;                    % free variables (to be populated)
    K(i).l = 0;                    % nonnegative variables (to be populated)
    K(i).s = [];                   % semidefinite variables (to be populated)
    K(i).q = 0;                    % to be ignored
    K(i).e = 0;                    % to be ignored
    K(i).c = 0;                    % to be ignored
    K(i).r = 0;                    % to be ignored
    K(i).p = 0;                    % to be ignored
    K(i).m = 0;                    % to be ignored
    K(i).scomplex = [];            % to be ignored
    K(i).xcomplex = [];            % to be ignored
    K(i).sos = [];                 % to be ignored
    K(i).schur_funs = [];          % to be ignored
    K(i).schur_data = [];          % to be ignored
    K(i).schur_variables = [];     % to be ignored
    
    % Then, empty cells At and c for the constraints in the moment relaxation.
    % Only consider equalities, linear inequalities (1D LMIs) and LMIs
    % c - At*y \in K.f (i.e., ==0)
    % c - At*y \in K.l (i.e., >=0)
    % c - At*y \in K.s (i.e., is PSD)
    At(i).f = {};
    At(i).l = {};
    At(i).s = {};
    c(i).f = {};
    c(i).l = {};
    c(i).s = {};
    
    % Keep track of indices in the clique
    momentID{i} = [];
    gramMonomials{i} = [];
    
    % Finally, the moments in this clique
    var_id = cliqueData.cliques{i}(:);
    M = monolistcoeff(length(var_id), 2*omega, 2*omega);
    [ii,jj,vv] = find(M(sum(M,2)>0, :));
    rows = [rows; ii+numRows];
    cols = [cols; var_id(jj)];
    vals = [vals; vv];
    numRows = numRows + max(ii);
end
all_moments = sparse(rows, cols, vals, numRows, numx);