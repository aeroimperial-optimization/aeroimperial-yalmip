function [r, S, U] = svd_rank(A, droptol)
% Find the rank of a matrix using drop in singular values as the deciding
% criterion. If no drop, we use the magnitude of the singular vales
% (1e-8 is used as a tolerance)
if nargout < 3
    S = svd(A);
    S = sort(S,'descend');
else
    [U,S] = svd(A,'econ','vector');
end

if all(S < 1e-8)
    r = 0;
else
    drop = S(2:end)./( eps + S(1:end-1) ); % avoid dividing by zero
    drop = find(drop < droptol,1,'first');
    if ~isempty(drop)
        r = drop;
    else
        r = nnz(S/S(1) >= 1e-8);
    end
end
end