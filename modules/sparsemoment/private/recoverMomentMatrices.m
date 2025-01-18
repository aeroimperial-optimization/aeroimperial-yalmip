function S = recoverMomentMatrices(sol)
%function S = recoverMomentMatrices(y, At, c, K, isMomentMatrix)

% Construct the moment matrices given the solution of a moment-SOS
% relaxation

% Empty input? Problem not solved, so return nothing
if isempty(sol.reducedMoments)
    S = [];
    return
end

% Slack variables and shift for constraints that comes before LMIs
s = sol.prog.c-sol.prog.At*sol.reducedMoments;
shift = sol.prog.K.f + sol.prog.K.l + sol.prog.K.q + sol.prog.K.r;

% Build the moment matrices
m = nnz(sol.prog.isMomentMatrix);
S = cell(m,1);
cnt = 1;
for i = 1:length(sol.prog.K.s)
    nsdp = sol.prog.K.s(i)^2;
    if sol.prog.isMomentMatrix(i)
        idx = shift + (1:nsdp);
        S{cnt} = reshape(s(idx), sol.prog.K.s(i), sol.prog.K.s(i));
        cnt  = cnt + 1;
    end
    shift = shift + nsdp;
end