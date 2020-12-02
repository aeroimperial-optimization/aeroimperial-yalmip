function cliques = corrSparsityCliques(x, OBJ, CNSTR, ordering)

% cliques = corrSparsityCliques(x, OBJ, CNSTR) Finds the cliques of the
%           chordal extension of the correlative sparsity graph for the 
%           given polynomial optimization problem (POP). Inputs:
%
%               x: the independent variables (sdpvar vector)
%               OBJ: the objective function of the POP (sdpvar)
%               CNSTR: the constraints of the POP (sdpvar vector)
%
%
% cliques = corrSparsityCliques(x, OBJ, CNSTR, ordering) also specifies
%           which matrix reordering is used to compute a chordal extension:
%
%               ordering = 0: symamd (default)
%               ordering = 1: symrcm

if nargin < 4
    ordering = 0;
else
   assert(ordering==0||ordering==1, ...
       'Error: type of ordering not supported (0 or 1)');
end

% Initialize
x = x(:);
n = numel(x);
C = eye(n,n);

% Get correlative sparsity matrix of objective:
% C(i,j) = 1 if variables x(i) and x(j) appear in a monomial
exponents = getexponentbase(OBJ, x);
for i = 1:size(exponents, 1)
    [~,var_idx] = find(exponents(i,:));
    for j = 1:length(var_idx)
        C(var_idx(j), var_idx(j)) = 1;
        for k = 2:length(var_idx)
            C(var_idx(j), var_idx(k)) = 1;
            C(var_idx(k), var_idx(j)) = 1;
        end
    end
end

% Get correlative sparsity matrix of constraints:
% C(i,j) = 1 if variables x(i) and x(j) appear in the same constraint.
% E.g. CNSTR = x1^2 + x1*x3 + x2 makes C(1:3,1:3) = 1.
for i = 1:length(CNSTR)
    exponents = getexponentbase(CNSTR(i), x);
    var_in_cnstr = sum(exponents, 1)>0;
    C(var_in_cnstr, var_in_cnstr) = 1;
end

% Make sure C is positive definite: ensure strict diagonal dominance
% Use 1-norm since matlab uses column-major storage.
C = C + norm(C,1).*eye(n,n);

% ----------------------------------------------------------------------- %
% Taken care of by cliquesFromSpMatD()?
% ----------------------------------------------------------------------- %
% % Chordal extension of sparsity pattern (use heuristic based on Cholesky)
% % Shift until cholesky does not fail
% s = symamd(C);
% shift = 100; p = 1;
% while p ~= 0
%     [R, p] = chol( C(s,s) + shift*speye(n) );
%     shift = 2*shift;
% end
% ----------------------------------------------------------------------- %

% Get cliques -- code from SparseCoLO
cliques = cliquesFromSpMatD(C,ordering);
cliques = cliques.Set(:);