function [sol] = admmSplitting(prog,options,initSol)

% [sol] = admmSplitting(prog,options)
% 
% Solve SDP relaxation of sparse POP using an ADMM method with clique
% decomposition strategy.
%
% Giovanni Fantuzzi
% Imperial College London
% April 2021

%============================================
% Import utils amd make options
tstart = tic;
import admm_utils.*
opts = admmOpts;
if(nargin >= 2); opts = setUserOpts(opts,options); end
if(nargin < 3); initSol = []; end

%============================================
% Print nice welcoming header
if opts.verbose
    [header,myline1,myline2] = printHeader();
    fprintf(myline1)
    fprintf('ADMM Solver with splitting (G. Fantuzzi)\n')
    fprintf(myline1)
    fprintf('Initializing solver...')
end

%============================================
% Initialize solver
proctime = tic;
[prog, opts] = checkInputs(prog, opts);
[sol, opts] = makeVariables(prog, opts, initSol);
[prog, opts, sol] = rescaleData(prog, opts, sol, ~isempty(initSol));
sol.time.setupTime = toc(proctime);
if opts.verbose
    fprintf('done in %.4f seconds.      \n',sol.time.setupTime);
    fprintf('Adaptive penalty       : %i\n',opts.adaptive);
    fprintf('Rescale data           : %i\n',opts.rescale);
    fprintf(myline1);
    fprintf(header);
    fprintf(myline2);
end

%============================================
% Run ADMM
sol.time.ADMMtime = tic;
iter = 0;
opts.stop = 0;
[sol, opts] = checkConvergence(prog, sol, opts, iter);
while ~opts.stop && iter <= opts.maxIter
    iter = iter + 1;
    sol = updateX(prog, sol, opts, iter);
    sol = updateY(prog, sol, opts, iter);
    sol = updateZ(prog, sol, opts, iter);
    [sol, opts] = checkConvergence(prog, sol, opts, iter);   
end
sol.time.totalTime = toc(tstart);
sol.penalty = opts.rho;

%============================================
% Print summary
if opts.verbose
    fprintf(myline1)
    fprintf(' SOLUTION SUMMARY:\n')
    fprintf('------------------\n')
    fprintf(' Number of iterations : %11.d\n',iter)
    fprintf(' Cost                 : %11.4e\n',sol.cost)
    fprintf(' Primal residual      : %11.4e\n',sol.pres)
    fprintf(' Dual residual        : %11.4e\n',sol.dres)
    fprintf('  Setup time   (s)    : %11.4e\n',sol.time.setupTime)
    fprintf('  ADMM  time   (s)    : %11.4e\n',sol.time.ADMMtime)
    fprintf(' Total  time   (s)    : %11.4e\n',sol.time.totalTime)
    fprintf('  Avg. block 1 (s)    : %11.4e\n',sol.time.updateX/iter)
    fprintf('  Avg. block 2 (s)    : %11.4e\n',sol.time.updateY/iter)
    fprintf(myline1)
end
