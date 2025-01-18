function x = extractminimizers(sol, droptol)

%function x = extractminimizers(momentMatrices, gramMonomials, droptol)

% Attempt to extract minimizers from a set of moment matrices. Inputs are:
% - momentMatrices: cell array of moment matrices
% -  gramMonomials: cell array of the corresponding monomial basis for the
%                   Gram matrix decomposition
%
% NOTE: We should check if the flatness condition holds, but we can just as
% well attempt to extract the minimizers directly.

if nargin < 2; droptol = 1e-3; end
momentMatrices = recoverMomentMatrices(sol);
if ~iscell(momentMatrices)
    momentMatrices = {momentMatrices};
    sol.prog.gramMonomials = {sol.prog.gramMonomials};
end

% Loop over moment matrices
for k = 1:length(momentMatrices)
    % Rank via SVD decomposition
    [rankM, S, U] = svd_rank(momentMatrices{k}, droptol);
    U = U(:,1:rankM)*diag(sqrt(S(1:rankM)));

    % Get column echelon form, like gloptipoly
    [U,basis] = cef(U,1e-6);

    % Decompose the powers in v and find monomials of degree 1
    beta = sol.prog.gramMonomials{k}(basis,:);
    nmons = size(sol.prog.gramMonomials{k},2);
    N = zeros(rankM^2,nmons);
    N2 = zeros(rankM^2,nmons);
    for j = 1:nmons
        pow = beta;
        pow(:,j) = pow(:,j) + 1;
        [ind,pos] = ismember(sol.prog.gramMonomials{k},pow,'rows');
        ind = find(ind);
        pos = pos(ind);
        T = zeros(rankM);
        T(pos,:) = U(ind,:);
        N(:,j) = T(:);
    end

    % Random convex combination of the (vectorized) matrices N(:,i)
    w = rand(nmons,1);
    w = w./sum(w);
    Y = reshape(N*w,rankM,rankM);

    % Order shur decomposition and extract minimizers
    [Q,T] = schur(Y);

    if isempty(T)
        % complex eigenvalues = failure
        x = [];

    else
        % Retrieve optimal vectors
        % It is assumed than there is no multiple root
        x{k} = zeros(nmons,rankM);
        for i = 1:rankM
            for j = 1:nmons
                x{k}(j,i) = Q(:,i)'*reshape(N(:,j),rankM,rankM)*Q(:,i);
            end
        end
    end

end

% else
%     % Make a cell and call again
%     x = extractminimizers({momentMatrices},{gramMonomials});
% end







