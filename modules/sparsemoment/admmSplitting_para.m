function [ws] = admmSplitting_para(prog,options,initSol)

% [sol] = admmSplitting_para(prog,options)
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
import admm_utils_para.*
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
[ws, opts] = checkInputs(prog, opts);
[ws, opts] = makeVariables(ws, opts, initSol);
ws(end).time.setupTime = toc(proctime);
if opts.verbose
    fprintf('done in %.4f seconds.      \n',ws(end).time.setupTime);
    fprintf('Adaptive penalty       : %i\n',opts.adaptive);
    fprintf('Rescale data           : %i\n',opts.rescale);
    fprintf(myline1);
    fprintf(header);
    fprintf(myline2);
end

%============================================
% Run ADMM
ws(end).time.ADMMtime = tic;
iter = 0;
opts.stop = 0;
[ws, opts] = checkConvergence(ws, opts, iter);
while ~opts.stop && iter <= opts.maxIter
    iter = iter + 1;
    ws = updateX(ws, opts, iter);
    ws = updateY(ws, opts, iter);
    ws = updateZ(ws, opts, iter);
    [ws, opts] = checkConvergence(ws, opts, iter);   
end
ws(end).time.totalTime = toc(tstart);
ws(end).penalty = opts.rho;

%============================================
% Print summary
if opts.verbose
    fprintf(myline1)
    fprintf(' SOLUTION SUMMARY:\n')
    fprintf('------------------\n')
    fprintf(' Number of iterations : %11.d\n',iter)
    fprintf(' Cost                 : %11.4e\n',ws(end).cost)
    fprintf(' Primal residual      : %11.4e\n',ws(end).pres)
    fprintf(' Dual residual        : %11.4e\n',ws(end).dres)
    fprintf('  Setup time   (s)    : %11.4e\n',ws(end).time.setupTime)
    fprintf('  ADMM  time   (s)    : %11.4e\n',ws(end).time.ADMMtime)
    fprintf(' Total  time   (s)    : %11.4e\n',ws(end).time.totalTime)
    fprintf('  Avg. block 1 (s)    : %11.4e\n',ws(end).time.updateX/iter)
    fprintf('  Avg. block 2 (s)    : %11.4e\n',ws(end).time.updateY/iter)
    fprintf(myline1)
end
