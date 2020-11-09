function [At, b, c, K,  unique_exponents] = makeConicUniqueMomentsNew(AA, b, c, K, exponents, shift, whichClique)

% Get unique exponents and combine cells in AA to build the constraint
% matrix At. Use sparse indexing for speed.
% Match cliques to powers of each clique in "exponent{i+2}"
% NOTE:
% exponents{1} = exponents in objective
% exponents{2} = exponents in linear constraints (all zeros)
% exponents{i+2} = exponents in clique i

% Initalize
% counter_jj = shift;
counter_ii = 0;
counter_vv = 0;
rows.f = cellfun(@(X)size(X, 1), AA.f);
rows.l = cellfun(@(X)size(X, 1), AA.l);
rows.s = cellfun(@(X)size(X, 1), AA.s);
all_rows = [rows.f, rows.l, rows.s];
cols.f = cellfun(@(X)size(X, 2), AA.f);
cols.l = cellfun(@(X)size(X, 2), AA.l);
cols.s = cellfun(@(X)size(X, 2), AA.s);
nonzeros.f = cellfun(@nnz, AA.f);
nonzeros.l = cellfun(@nnz, AA.l);
nonzeros.s = cellfun(@nnz, AA.s);
kk = sum(nonzeros.f) + sum(nonzeros.l) + sum(nonzeros.s);
II = zeros(kk,1);
JJ = zeros(kk,1);
VV = zeros(kk,1);

% Gen # moments for each clique
num_mons_clique = cellfun(@(X)size(X,1),exponents);
total_shifts = cumsum([1, num_mons_clique])-1;

% Get unique exponents (overwrite exponents)
unique_exponents = vertcat(exponents{:});
[~, IA, IC] = unique(unique_exponents*rand(size(unique_exponents,2),1));
unique_exponents = unique_exponents(IA,:);
m = size(unique_exponents, 1);

% Make objective
b = accumarray(IC(1:shift), b, [m, 1]);

% Assemble At
num_free = length(AA.f);
num_lin  = length(AA.l);
num_sdp  = length(AA.s);
num_cnstr = num_free + num_lin + num_sdp;
for k = 1:num_cnstr
    if k < num_free
        % free cone!
        idx = k;
        clique_idx = whichClique.f(idx)+2;
        idx_jj = total_shifts(clique_idx) + (1:cols.f(idx));
        idx_vv = counter_vv + (1:nonzeros.f(idx));
        localIC = IC(idx_jj);
        [iA, jA, VV(idx_vv)] = find(AA.f{idx});
        II(idx_vv) = iA + counter_ii;
        JJ(idx_vv) = localIC(jA);
        counter_ii = counter_ii + rows.f(idx);
        counter_vv = counter_vv + nonzeros.f(idx);
        
    elseif k == num_free
        % free cone -- mass constraint
        clique_idx = 2;
        idx = k;
        idx_jj = total_shifts(clique_idx) + (1:cols.f(idx));
        idx_vv = counter_vv + (1:nonzeros.f(idx));
        localIC = IC(idx_jj);
        [iA, jA, VV(idx_vv)] = find(AA.f{idx});
        II(idx_vv) = iA + counter_ii;
        JJ(idx_vv) = localIC(jA);
        counter_ii = counter_ii + rows.f(idx);
        counter_vv = counter_vv + nonzeros.f(idx);
        
    elseif (k > num_free) && (k <= num_free + num_lin)
        % linear cone
        % Nothing to do!
        
    elseif (k > num_free + num_lin) && (k <= num_free + num_lin + num_sdp)
        % SDP cone
        idx = k - num_free - num_lin;
        clique_idx = whichClique.s(idx)+2;
        idx_jj = total_shifts(clique_idx) + (1:cols.s(idx));
        idx_vv = counter_vv + (1:nonzeros.s(idx));
        localIC = IC(idx_jj);
        [iA, jA, VV(idx_vv)] = find(AA.s{idx});
        II(idx_vv) = iA + counter_ii;
        JJ(idx_vv) = localIC(jA);
        counter_ii = counter_ii + rows.s(idx);
        counter_vv = counter_vv + nonzeros.s(idx);
        
    end
    
end
At = sparse(II,JJ,VV,sum(all_rows),m);

% Make c: simple concatenation of vectors
c = [c.f; c.l; c.s];

% End function
end
