function [At, b, c, K,  exponents] = makeConicUniqueMoments(AA, b, c, K, exponents, shift)

% Get unique exponents and combine cells in AA to build the constraint
% matrix At. Use sparse indexing for speed.

% Initalize
counter_jj = shift;
counter_ii = 0;
counter_vv = 0;
rows = cellfun(@(X)size(X, 1), AA);
cols = cellfun(@(X)size(X, 2), AA);
nonzeros = cellfun(@nnz, AA);
kk = sum(nonzeros);
II = zeros(kk,1);
JJ = zeros(kk,1);
VV = zeros(kk,1);

% Get unique exponents
% [exponents, ~, IC] = unique(exponents, 'rows');
[~, IA, IC] = unique(exponents*rand(size(exponents,2),1));
exponents = exponents(IA,:);
m = size(exponents, 1);
b = accumarray(IC, b, [m, 1]);

% Assemble At
for k = 1:length(AA)
    idx_jj = counter_jj + (1:cols(k));
    idx_vv = counter_vv + (1:nonzeros(k));
    localIC = IC(idx_jj);
    [iA, jA, VV(idx_vv)] = find(AA{k});
    II(idx_vv) = iA + counter_ii;
    JJ(idx_vv) = localIC(jA);
    counter_jj = idx_jj(end);
    counter_ii = counter_ii + rows(k);
    counter_vv = counter_vv + nonzeros(k);
end

At = sparse(II,JJ,VV,sum(rows),m);

% End function
end
