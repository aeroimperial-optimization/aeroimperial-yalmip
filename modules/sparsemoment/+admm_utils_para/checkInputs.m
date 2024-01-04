function [ws, opts] = checkInputs(prog,opts)

% CHECKINPUTS
% Verify that inputs are consistent, extract some useful variables, and
% create the workspace

% Correct sizes?
opts.noCliques = length(prog.At);
opts.m.y = length(prog.b.y);
assert(length(prog.b.s)==opts.noCliques, 'Length of input (cell) arrays must match')
assert(length(prog.c)==opts.noCliques, 'Length of input (cell) arrays must match')
assert(length(prog.d)==opts.noCliques, 'Length of input (cell) arrays must match')
assert(length(prog.E)==opts.noCliques, 'Length of input (cell) arrays must match')
assert(length(prog.F)==opts.noCliques, 'Length of input (cell) arrays must match')
assert(length(prog.K)==opts.noCliques, 'Length of input (cell) arrays must match')

% Set up the workspace: an array with slices corresponding to cliques
% Cone is copied later
[ws(1:opts.noCliques).At] = deal(prog.At{:});
[ws(1:opts.noCliques).bs] = deal(prog.b.s{:});
[ws(1:opts.noCliques).c] = deal(prog.c{:});
[ws(1:opts.noCliques).d] = deal(prog.d{:});
[ws(1:opts.noCliques).E] = deal(prog.E{:});
[ws(1:opts.noCliques).F] = deal(prog.F{:});

% Loop over cliques, assign cone, and check inter-clique dimensions
opts.m.s = zeros(opts.noCliques,1);
opts.n.s = zeros(opts.noCliques,1);
for i = 1:opts.noCliques
    [opts.n.s(i), opts.m.s(i)] =  size(ws(i).At);
    [ws(i).K, nConeVars] = setCone(prog.K(i));
    assert(length(ws(i).bs)==opts.m.s(i), 'Data dimensions mismatch')
    assert(length(ws(i).c)==opts.n.s(i), 'Data dimensions mismatch')
    assert(nConeVars==opts.n.s(i), 'Data dimensions mismatch')
    assert(size(ws(i).E,2)==opts.m.s(i), 'Data dimensions mismatch')
    assert(size(ws(i).F,2)==opts.m.y, 'Data dimensions mismatch')
    assert(size(ws(i).F,1)==size(ws(i).E,1), 'Data dimensions mismatch')
    assert(size(ws(i).F,1)==length(ws(i).d), 'Data dimensions mismatch')
    assert(size(ws(i).E,1)==length(ws(i).d), 'Data dimensions mismatch')
end

% Make last component of the workspace: vector y
ws(opts.noCliques+1).b = prog.b.y;

% Compute "complement" of lambda (used later to compute cost) and scaled
% parameters (used in ADMM iterations)
opts.lambdaComp = 1-opts.lambda;
opts.lambdaOverRho = opts.lambda / opts.rho;
opts.lambdaCompOverRho = opts.lambdaComp / opts.rho;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, nConeVars] = setCone(K)
% Initialize cone
nConeVars = 0;
% Free variables?
if(isfield(K,'f') && ~isempty(K.f) && K.f > 0)
    nConeVars = nConeVars + K.f;
else
    K.f = 0;
end
% Nonnegative variables?
if(isfield(K,'l') && ~isempty(K.l) && K.l > 0)
    nConeVars = nConeVars + K.l;
else
    K.l = 0;
end
% Quadratic cones?
if (isfield(K,'q') && ~isempty(K.q) && max(K.q) > 0)
    K.q = K.q(K.q~=0);
    nConeVars = nConeVars + sum(K.q);
else
    K.q = 0;
end
% PSD cones?
if (isfield(K,'s') && ~isempty(K.s) && max(K.s) > 0)
    K.s = K.s(K.s~=0);
    nConeVars = nConeVars + sum(K.s.^2);
else
    K.s = 0;
end
end