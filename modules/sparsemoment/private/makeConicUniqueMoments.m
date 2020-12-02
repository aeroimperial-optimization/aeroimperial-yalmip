function [At, b_in, c, K] = makeConicUniqueMoments(At_in, b_in, c_in, K_in)

% Combine constraints on individual cliques into total conic program.
% Inputs At, c and K are structures with conic constraints for each clique.
% In this function, b is already assembled so there is no need to operate
% on it.

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
temp = [At_in(:).f]; At = vertcat(temp{:});
temp = [c_in(:).f]; c = vertcat(temp{:});
K.f = size(At,1);

% Assemble the linear constraints
temp = [At_in(:).l]; temp = vertcat(temp{:}); At = vertcat(At, temp);
temp = [c_in(:).l]; temp = vertcat(temp{:}); c = vertcat(c, temp);
K.l = size(temp,1);

% Assemble the semidefinite constraints
temp = [At_in(:).s]; At = vertcat(At, temp{:});
temp = [c_in(:).s]; c = vertcat(c, temp{:});
K.s = [K_in(:).s];

% End function
end
