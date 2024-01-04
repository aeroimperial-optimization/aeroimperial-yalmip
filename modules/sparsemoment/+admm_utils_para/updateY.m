function ws = updateY(ws, opts, iter)

% UPDATEY
% This function updates sol.y and sol.z{i} for all indices i
% Solution of y is cheap: diagonal system (but assempling RHS can be slow?)
% Solution of z{i} requires a conic projection
% At first iteration, we set up useful persistent data

% Initialize
tstart = tic;
if iter==1
    P = zeros(opts.m.y,1);
    for i = 1:opts.noCliques
        P = P + diag(ws(i).F.' * ws(i).F);
        ws(i).dresvec = zeros(size(ws(i).At, 2), 1);
        ws(i).dres = 0;
    end
    ws(end).P = 1./P;
end

% Update y and z(i) in a single loop
RHS = opts.lambdaOverRho .* ws(end).b;
parfor i = 1:opts.noCliques
    RHS = RHS + ws(i).F.' * ( ws(i).d + ws(i).xi - ws(i).E*ws(i).s );
    v = ws(i).c - ws(i).At * ws(i).s + ws(i).eta;
    zNew = projectK(v, ws(i).K);
    ws(i).dresvec(:) = ws(i).At.' * (zNew - ws(i).z);
    ws(i).z = zNew;
end
yNew = RHS .* ws(end).P;
yDiff = yNew - ws(end).y;
ws(end).y = yNew;

% Update dual residual (another loop needed? For now. Can do better by splitting the y update across iterations)
parfor i = 1:opts.noCliques
    temp = ws(i).F * yDiff;
    temp = ws(i).E.' * temp;
    ws(i).dresvec(:) = ws(i).dresvec(:) +  temp;
    ws(i).dres = norm(ws(i).dresvec, 'inf');
end

% Finished
ws(end).dres = opts.rho .* max([ws(1:opts.noCliques).dres]);
ws(end).time.updateY = ws(end).time.updateY + toc(tstart);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = projectK(X,K)
shift = 0;
% ZERO CONE PROJECTION
if(isfield(K,'f') && K.f > 0)
    pos = shift + (1:K.f);
    X(pos) = projectZero(X(pos));
    shift = shift + K.f;
end
% POSITIVE ORTHANT PROJECTION
if(isfield(K,'l') && K.l > 0)
    pos = shift + (1:K.l);
    X(pos) = projectPosOrthant(X(pos));
    shift = shift + K.l;
end
% SOC PROJECTION
if (isfield(K,'q') && any(K.q))
    for i = (1:length(K.q))
        pos = shift + (1:K.q(i));
        X(pos) = projectSOC(X(pos));
        shift = shift + K.q(i);
    end
end
% PSD CONE PROJECTION
if (isfield(K,'s') && any(K.s))
    for i = (1:length(K.s))
        pos = shift + (1:K.s(i)^2);
        X(pos) = projectPSD(X(pos), K.s(i));  %self dual cone
        shift = shift + K.s(i)^2;
    end
end
end

%------------------------------------
% ACTUAL PROJECTION FUNCTIONS
%------------------------------------
% ZERO CONE
function X = projectZero(X)
X(:) = 0;
end

% POSITIVE ORTHANT
function X = projectPosOrthant(X)
X = max(X,0);
end

% SECOND ORDER CONE
function X = projectSOC(X)
s = X(1);
x = X(2:end);
normx = norm(x,2);
if(normx <= s)
    X = X;
elseif(normx <= -s)
    X = X.*0;
else
    t = (normx + s)/2;
    X = [t;(t/normx).*x];
end
end

% POSITIVE SEMIDEFINITE CONE
function s = projectPSD(s, n)
S = blockify(zeros(n), 0.5.*s);
S = S+S.';
[U,E] = eig(S);
ind = diag(E)>0;
UsE = U(:,ind)*sqrt(E(ind,ind));
S = UsE*UsE.';
s(:) = S(:);
end

% VECTOR TO MATRIX
function X = blockify(X,x)
% Reshape vector obtained with vec(X) into the matrix X
X(:) = x;
end