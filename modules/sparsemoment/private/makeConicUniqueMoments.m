function [At, b_in, c, K, isMomentMatrix, PROJ, bshift, y0] = makeConicUniqueMoments(At_in, b_in, c_in, K_in, isMomentMatrix, options)

% Combine constraints on individual cliques into total conic program.
% Inputs At, c and K are structures with conic constraints for each clique.
% Also try to eliminate known moments and represent original moment vector
% as y = y0+ PROJ*z

% First, an empty cone with all fields required by YALMIP
K.f = 0;                    % free variables (to be populated)
K.l = 0;                    % nonnegative variables (to be populated)
K.s = [];                   % semidefinite variables (to be populated)
K.q = 0;                    % to be ignored
K.e = 0;                    % to be ignored
K.c = 0;                    % to be ignored
K.r = 0;                    % to be ignored
K.p = 0;                    % to be ignored
K.m = 0;                    % to be ignored
K.scomplex = [];            % to be ignored
K.xcomplex = [];            % to be ignored
K.sos = [];                 % to be ignored
K.schur_funs = [];          % to be ignored
K.schur_data = [];          % to be ignored
K.schur_variables = [];     % to be ignored

% Assemble the free cone
temp = [At_in(:).f]; At_free = vertcat(temp{:});
temp = [c_in(:).f]; c_free = vertcat(temp{:});
K.f = size(c_free,1);

% Assemble the linear constraints
temp = [At_in(:).l]; At_lin = vertcat(temp{:});
temp = [c_in(:).l]; c_lin = vertcat(temp{:});
K.l = size(c_lin,1);

% Assemble the semidefinite constraints
temp = [At_in(:).s]; At_sdp = vertcat(temp{:});
temp = [c_in(:).s]; c_sdp = vertcat(temp{:});
K.s = [K_in(:).s];

% Do we have any 1-by-1 PSD cones? Move to linear cones
% We assume moment matrices are not 1-by-1...
islin = (K.s==1);
if any(islin)
    K.l = K.l + sum(islin);
    idx = cumsum(K.s.^2);
    K.s = K.s(~islin);
    isMomentMatrix = isMomentMatrix(~islin);
    rows = 1:idx(end);
    islin = ismember(rows, idx(islin));
    At_lin = [At_lin; At_sdp(islin,:)];
    c_lin = [c_lin; c_sdp(islin)];
    At_sdp = At_sdp(~islin,:);
    c_sdp = c_sdp(~islin);
end

% Assemble
At = [At_free; At_lin; At_sdp];
c = [c_free; c_lin; c_sdp];

% Remove simple equalities:
% (1) moments with fixed values
% (2) moments that are proportional to one other moment
% And iterate to remove as many moments as possible. The original moments
% are recovered by y = y0 + PROJ*z and the original cost is
if options.sparsemoment.eliminateMoments
    [At,b_in,c,K,PROJ,bshift,y0] = eliminateMatchingMoments(At,b_in,c,K);
else
    num_moments = size(At,2);
    PROJ = speye(num_moments);
    bshift = 0;
    y0 = sparse(num_moments,1);
end


% End function
end
