function prog = makeConicUniqueMoments(At_in, b_in, c_in, K_in, isMomentMatrix, options)

% Combine constraints on individual cliques into total conic program.
% Inputs At, c and K are structures with conic constraints for each clique.
% Also try to eliminate known moments and represent original moment vector
% as y = y0+ PROJ*z

% Initialize program
prog.At = [];
prog.b = b_in;
prog.c = [];
prog.K = [];
prog.isMomentMatrix = horzcat(isMomentMatrix{:});
prog.PROJ = [];
prog.bshift = [];
prog.y0 = [];

% First, an empty cone with all fields required by YALMIP
prog.K.f = 0;                    % free variables (to be populated)
prog.K.l = 0;                    % nonnegative variables (to be populated)
prog.K.s = [];                   % semidefinite variables (to be populated)
prog.K.q = 0;                    % to be ignored
prog.K.e = 0;                    % to be ignored
prog.K.c = 0;                    % to be ignored
prog.K.r = 0;                    % to be ignored
prog.K.p = 0;                    % to be ignored
prog.K.m = 0;                    % to be ignored
prog.K.scomplex = [];            % to be ignored
prog.K.xcomplex = [];            % to be ignored
prog.K.sos = [];                 % to be ignored
prog.K.schur_funs = [];          % to be ignored
prog.K.schur_data = [];          % to be ignored
prog.K.schur_variables = [];     % to be ignored

% Assemble the free cone
temp = [At_in(:).f]; At_free = vertcat(temp{:});
temp = [c_in(:).f]; c_free = vertcat(temp{:});
prog.K.f = size(c_free,1);

% Assemble the linear constraints
temp = [At_in(:).l]; At_lin = vertcat(temp{:});
temp = [c_in(:).l]; c_lin = vertcat(temp{:});
prog.K.l = size(c_lin,1);

% Assemble the semidefinite constraints
temp = [At_in(:).s]; At_sdp = vertcat(temp{:});
temp = [c_in(:).s]; c_sdp = vertcat(temp{:});
prog.K.s = [K_in(:).s];

% Do we have any 1-by-1 PSD cones? Move to linear cones
% We assume moment matrices are not 1-by-1...
islin = (prog.K.s==1);
if any(islin)
    prog.K.l = prog.K.l + sum(islin);
    idx = cumsum(prog.K.s.^2);
    prog.K.s = prog.K.s(~islin);
    prog.isMomentMatrix = prog.isMomentMatrix(~islin);
    rows = 1:idx(end);
    islin = ismember(rows, idx(islin));
    At_lin = [At_lin; At_sdp(islin,:)];
    c_lin = [c_lin; c_sdp(islin)];
    At_sdp = At_sdp(~islin,:);
    c_sdp = c_sdp(~islin);
end

% Assemble
prog.At = [At_free; At_lin; At_sdp];
prog.c = [c_free; c_lin; c_sdp];

% Remove simple equalities:
% (1) moments with fixed values
% (2) moments that are proportional to one other moment
% And iterate to remove as many moments as possible. The original moments
% are recovered by y = y0 + PROJ*z and the original cost is
if options.sparsemoment.eliminateMoments
    [prog.At,prog.b,prog.c,prog.K,prog.PROJ,prog.bshift,prog.y0] = eliminateMatchingMoments(prog.At,prog.b,prog.c,prog.K);
else
    num_moments = size(prog.At,2);
    prog.PROJ = speye(num_moments);
    prog.bshift = 0;
    prog.y0 = sparse(num_moments,1);
end


% End function
end
