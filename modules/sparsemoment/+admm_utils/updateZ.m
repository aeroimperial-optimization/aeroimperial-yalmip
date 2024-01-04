function sol = updateZ(prog, sol, opts, iter)

% UPDATEZ
% Update the multipliers

% Initalize some useful variables to update the residual
tstart = tic;
persistent pnorm
if iter==1 || isempty(pnorm)
   pnorm = zeros(2*opts.noCliques,1);
end

% Operate
for i = 1:opts.noCliques
    p1 = prog.c{i} - prog.At{i} * sol.s{i} - sol.z{i};
    p2 = prog.d{i} - prog.E{i} * sol.s{i} - prog.F{i} * sol.y;
    sol.eta{i} = sol.eta{i} + p1;
    sol.xi{i} = sol.xi{i} + p2;
    pnorm(i) = norm(p1,'Inf');
    pnorm(opts.noCliques+i) = norm(p2,'Inf');
end

% Update the primal residual
sol.pres = max(pnorm);
sol.time.updateZ = sol.time.updateZ + toc(tstart);