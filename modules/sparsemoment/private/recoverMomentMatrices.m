function S = recoverMomentMatrices(y, At, c, K, isMomentMatrix)

% Construct the moment matrices given the solution of a moment-SOS
% relaxation

% Slack variables and shift for constraints that comes before LMIs
s = c-At*y;
shift = K.f + K.l + K.q + K.r;

% Build the moment matrices
m = nnz(isMomentMatrix);
S = cell(m,1);
cnt = 1;
for i = 1:length(K.s)
    nsdp = K.s(i)^2;
    if isMomentMatrix(i)
        idx = shift + (1:nsdp);
        S{cnt} = reshape(s(idx), K.s(i), K.s(i));
        cnt  = cnt + 1;
    end
    shift = shift + nsdp;
end