function [prog, opts, sol] = rescaleData(prog, opts, sol, usrInitSol)

% +admm_utils/RESCALEDATA.m
% Try to rescale data to get nicer convergence properties.
% The equality constraints have the form
%
% [c{1}]   [At{1}   0    .....    0   | I  0 ... 0 |  0  ]   [s{1}]   [0]
% [d{1}]   [ E{1}   0    .....    0   | 0  0 ... 0 | F{1}]   [....]   [0]
% [c{2}]   [  0   At{2}  .....    0   | I  0 ... 0 |  0  ]   [s{k}]   [0]
% [d{2}] - [  0    E{2}  .....    0   | 0  0 ... 0 | F{2}] * [z{1}] = [0]
% [....]   [                  ....                       ]   [....]   [0]
% [c{k}]   [..... .....  .....  At{k} | 0  0 ... I |  0  ]   [z{k}]   [0]
% [d{k}]   [..... .....  .....   E{k} | 0  0 ... 0 | F{k}]   [ y  ]   [0]
%
% Using the same ordering of the variables, the cost vector is
%
% [b.s{1} ; ... ; b.s{k} ; 0 ; ... ; 0 ; b.y]
%
% We want:
% (1) Columns of the matrix to have norm similar to the norm of the constant vector
% (2) Rows of the matrix and b should have similar norm

% Parameters
minmax_scales = [1, +Inf];
nPasses = 1;

% Initialize
opts.scaleFactors.y = 1; % same for all cliques!
[opts.scaleFactors.s{1:opts.noCliques}] = deal(1);
[opts.scaleFactors.eta{1:opts.noCliques}] = deal(1); % This is also the scaling factor for z{i}
[opts.scaleFactors.xi{1:opts.noCliques}] = deal(1);

% Compute
if opts.rescale

    % --------------------------------------------------------------------------
    % Rescale multiple times in a loop
    for k = 1:nPasses

        % Loop over cliques
        for i = 1:opts.noCliques

            % Balance rows of At{i} in constraint c{i} - At{i}*s{i} \in K{i}
            % This is the scaling factor for eta{i} and z{i}
            [mA,nA] = size(prog.At{i});
            minmaxScale = minmax_scales .*sqrt(mA+nA);
            D1 = rowScalingMatrix(prog.At{i}, prog.K(i), minmaxScale);
            prog.At{i} = spdiags(D1,0,mA,mA) * prog.At{i};
            opts.scaleFactors.eta{i} = D1 .* opts.scaleFactors.eta{i};

            % Balance rows of E{i} and F{i} in constraint d{i} - E{i}*s{i} - F{i}*y = 0
            % This is the scaling factor for xi{i}
            [mE,nE] = size(prog.E{i});
            minmaxScale = minmax_scales .*sqrt(mE+nE);
            K.f = mE; K.l = 0; K.q = []; K.s = [];
            D2 = rowScalingMatrix([prog.E{i}, prog.F{i}], K, minmaxScale);
            prog.E{i} = spdiags(D2,0,mE,mE) * prog.E{i};
            prog.F{i} = spdiags(D2,0,mE,mE) * prog.F{i};
            opts.scaleFactors.xi{i} = D2 .* opts.scaleFactors.xi{i};

            % Balance cols of At{i} and E{i} (scaling must be the same)
            % This is the scaling factor for s{i}
            minmaxScale = minmax_scales .*sqrt(mA+mE+nE);
            E = colScalingMatrix([prog.At{i}; prog.E{i}], [], minmaxScale);
            prog.At{i} = prog.At{i} * spdiags(E,0,nE,nE);
            prog.E{i} = prog.E{i} * spdiags(E,0,nE,nE);
            opts.scaleFactors.s{i} = opts.scaleFactors.s{i} .* E;

            % Scale b.s{i}, c{i} and d{i}
            prog.b.s{i} = prog.b.s{i} .* E;
            prog.c{i} = prog.c{i} .* D1;
            prog.d{i} = prog.d{i} .* D2;
        end

        % Rescale cols of F{i}: must scale each i with the same value
        % This is the scaling factor for y
        mF = sum( cellfun(@(X)size(X,1), prog.F) );
        nF = size(prog.F{1}, 2);
        minmaxScale = minmax_scales .*sqrt(mF+nF);
        E = colScalingMatrix(vertcat(prog.F{i}), [], minmaxScale);
        opts.scaleFactors.y = opts.scaleFactors.y .* E;
        H = spdiags(E,0,nF,nF);
        for i = 1:opts.noCliques
            prog.F{i} = prog.F{i} * H;
        end

        % Rescale b.y
        prog.b.y = prog.b.y .* E;
    end

    % --------------------------------------------------------------------------
    % Now rescale the initial solution (account for user-specific solution)
    if usrInitSol
        sol.y = sol.y ./ opts.scaleFactors.y;
        for i = 1:opts.noCliques
            sol.s{i}   = sol.s{i}   ./ opts.scaleFactors.s{i};
            sol.z{i}   = blockify( sol.z{i}, admm_utils.flatten(sol.z{i}, prog.K(i))   ./ opts.scaleFactors.eta{i}, prog.K(i) );
            sol.eta{i} = sol.eta{i} ./ opts.scaleFactors.eta{i};
            sol.xi{i}  = sol.xi{i}  ./ opts.scaleFactors.xi{i};
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = colScalingMatrix(At, K, minmaxScale)
% Scaling of columns for the constraint c-At*s \in K
% Here K is not used because s is a free variables and its entries can be
% scaled arbitrarily.
E = full(sqrt(sum(At.^2,1)));  % norm of cols of A (row vec)
E(E>minmaxScale(2)) = minmaxScale(2);
E(E<minmaxScale(1)) = minmaxScale(1);
% Invert (taking possible zeros into account)
E = E(:);
ind = E~=0;
E(ind) = 1./E(ind);
E(~ind) = 1;
end

function D = rowScalingMatrix(At, K, minmaxScale)
% Construct matrix for row scaling in constraint c-At*s \in K
% Need to ensure that all rows corresponding to a single constraints are
% scaled by the same constant. This is important for quadratic and SDP
% cones only (the free/zero cone and the nonnegative orthant are fine)
D = full(sqrt(sum(At.^2,2)));  % norm of rows of A (col vec)
count = K.f + K.l; % Shift by free and linear cones
% Average of rows on quadratic cones
if sum(K.q)>0
    for i = 1:length(K.q)
        nvars = K.q(i);
        D(count+1:count+nvars) = mean(D(count+1:count+nvars));
        count = count + nvars;
    end
end
% Average of rows on PSD cones
if sum(K.s)>0
    for i = 1:length(K.s)
        nvars = K.s(i)^2;
        D(count+1:count+nvars) = mean(D(count+1:count+nvars));
        count = count + nvars;
    end
end
% Impose bounds on scaling
D(D>minmaxScale(2)) = minmaxScale(2);
D(D<minmaxScale(1)) = minmaxScale(1);
% Invert (taking possible zeros into account)
D = D(:);
ind = (D~=0);
D(ind) = 1./D(ind);
D(~ind) = 1;
end