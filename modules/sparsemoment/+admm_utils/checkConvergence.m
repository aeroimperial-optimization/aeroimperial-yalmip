function [sol, opts] = checkConvergence(prog, sol, opts, iter)

% CHECKCONVERGENCE
% Check termination criteria of ADMM method

% initialize persistent iteration counters when entering for the first
% time (persistent variables are empty and not zero when declared)

% Update cost
sol.cost = -opts.lambda * ( prog.b.y.' * sol.y );
for i = 1:opts.noCliques
    sol.cost = sol.cost - opts.lambdaComp * ( prog.b.s{i}.' * sol.s{i} );
end

%stopping criteria
runtime = toc(sol.time.ADMMtime);
if ( sol.pres<opts.relTol && sol.dres<opts.relTol ) || (iter==opts.maxIter)
    opts.stop = true;
    sol.time.ADMMtime = runtime;
end

%progress message
if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || opts.stop)
    fprintf('%5d | %8.2e | %8.2e | %9.2e  | %8.2e | %8.2e |\n',...
        iter,sol.pres,sol.dres,sol.cost,opts.rho,runtime);
end

% Update penalty parameter
if opts.adaptive; [sol, opts] = updatePenaltyParameter(iter, sol, opts); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sol, opts] = updatePenaltyParameter(iter, sol, opts)
% Update penalty parameter and take into account rescaling of the dual variables
persistent itPinf itDinf
if iter == 1 || isempty(itPinf)
    itPinf = 0; % number of iterations for which pinf/dinf <= eta
    itDinf = 0; % number of iterations for which pinf/dinf > eta
end
update = false;
resRat = sol.pres/sol.dres;
if resRat >= opts.nu
    itPinf = itPinf+1;
    itDinf = 0;
    if itPinf >= opts.rhoIt
        % ratio of pinf and dinf remained large for long => rescale rho
        itPinf = 0;
        newRho = min(opts.rho*opts.mu, opts.rhoMax);
        update = true;
    end
elseif 1/resRat >= opts.nu
    itDinf = itDinf+1;
    itPinf = 0;
    if itDinf >= opts.rhoIt
        % ratio of pinf and dinf remained small for long => rescale rho
        itDinf = 0;
        newRho = max(opts.rho/opts.mu, opts.rhoMin);
        update = true;
    end
end


if update
    scaling = opts.rho/newRho;
    % TODO: Check which one is faster!
%     sol.eta = cellfun(@(X)scaling.*X, sol.eta, 'UniformOutput', 0);
%     sol.xi = cellfun(@(X)scaling.*X, sol.xi, 'UniformOutput', 0);
    for i = 1:opts.noCliques
        sol.eta{i} = scaling .* sol.eta{i};
        sol.xi{i} = scaling .* sol.xi{i};
    end
    opts.rho = newRho;
    opts.lambdaComp = 1-opts.lambda;
    opts.lambdaOverRho = opts.lambda / opts.rho;
    opts.lambdaCompOverRho = opts.lambdaComp / opts.rho;
end
end