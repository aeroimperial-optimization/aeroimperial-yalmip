function [success, r, RR] = checkFlatExtension(sol, tol)

% Check if flat extension conditions hold for sparse POP with given
% cliques. The input sol is the solution structure returned by
% solvesparsemoment, while droptol is an optional tolerance for rank
% detection using an SVD-based criterion.

% The flat extension conditions are a refinement of the results given in:
% Theorem 
% Theorem 3.4 of  https://arxiv.org/abs/2406.06882

% Check arguments
if nargin < 2; tol = 1e-3; end

% Initalize
success = true; % innocent until proven guilty
r = NaN(sol.cliques.NoC, 1);
d = sol.relax_order - sol.relax_order_cnstr;

% rankfun = @(A) rank(A, 1e-8);
rankfun = @(A) svd_rank(A, tol);

% Recover the moment matrices
momentMatrices = recoverMomentMatrices(sol);

% Check flat extension conditions
RR = NaN(sol.cliques.NoC,sol.cliques.NoC);
try
    for k = 1 : sol.cliques.NoC
        % rank of moment matrix
        r(k) = rankfun(momentMatrices{k});
        % rank of submatrix of smaller relaxation order
        n_clique = numel(sol.cliques.Set{k});
        keep = sum(sol.prog.gramMonomials{k}, 2) <= d(k);
        r_new = rankfun(momentMatrices{k}(keep,keep));
        assert(r_new==r(k),'Rank Mismatch');
        % Check rank of submatrices of intersection between cliques
        % NOTE: It is enough to check for the cliques that we have not yet
        % looped over since the overlaps are identical.
        RR(k,k) = r(k);
        for j = k+1 : sol.cliques.NoC
            % Find moments supported on overlap
            [C,iA] = intersect(sol.cliques.Set{k}, sol.cliques.Set{j});
            if ~isempty(C)
                remove = true(n_clique,1); remove(iA) = false;
                overlap = sum(sol.prog.gramMonomials{k}(:,remove), 2) == 0;
                % Rank of overlap matrix
                r_overlap = rankfun(momentMatrices{k}(overlap,overlap));
                assert(r_overlap==r(k) || r_overlap==1,'Rank Mismatch');
                % Rank of overlap matrix of smaller degree
                keep = overlap & ( sum(sol.prog.gramMonomials{k}, 2) <= sol.relax_order-1 );
                r_new = rankfun(momentMatrices{k}(keep,keep));
                assert(r_new==r_overlap,'Rank Mismatch');
                RR(j,k) = r_overlap; RR(k,j) = RR(j,k);
            end
        end
    end
    assert(all(r==r(1)))
    r = r(1);
catch
    % We failed a rank test, return false
    success = false;
    r = [];
end

% END OF MAIN FUNCTION
end
