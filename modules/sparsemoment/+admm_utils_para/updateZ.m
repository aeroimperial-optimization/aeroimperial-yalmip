function ws = updateZ(ws, opts, iter)

% UPDATEZ
% Update the multipliers

% Initalize some useful variables to update the residual
tstart = tic;

% Operate
y = ws(end).y;
parfor i = 1:opts.noCliques
    p1 = ws(i).c - ws(i).At * ws(i).s - ws(i).z;
    p2 = ws(i).d - ws(i).E * ws(i).s - ws(i).F * y;
    ws(i).eta = ws(i).eta + p1;
    ws(i).xi = ws(i).xi + p2;
    ws(i).pres = norm(p1,'Inf');
    ws(i).pres = max(ws(i).pres, norm(p2,'Inf'));
end

% Update the primal residual
ws(end).pres = max([ws(1:opts.noCliques).pres]);
ws(end).time.updateZ = ws(end).time.updateZ + toc(tstart);