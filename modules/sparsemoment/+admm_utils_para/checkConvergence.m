function [ws, opts] = checkConvergence(ws, opts, iter)

% CHECKCONVERGENCE
% Check termination criteria of ADMM method

% initialize persistent iteration counters when entering for the first
% time (persistent variables are empty and not zero when declared)

% Update cost
ws(end).cost = -opts.lambda * ( ws(end).b.' * ws(end).y );
lambdaComp = opts.lambdaComp;
parfor i = 1:opts.noCliques
    ws(i).cost = -lambdaComp * ( ws(i).bs.' * ws(i).s );
end
ws(end).cost = sum([ws(:).cost]); 

%stopping criteria
runtime = toc(ws(end).time.ADMMtime);
if ( ws(end).pres<opts.relTol && ws(end).dres<opts.relTol ) || (iter==opts.maxIter)
    opts.stop = true;
    ws(end).time.ADMMtime = runtime;
end

%progress message
if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || opts.stop)
    fprintf('%5d | %8.2e | %8.2e | %9.2e  | %8.2e | %8.2e |\n',...
        iter,ws(end).pres,ws(end).dres,ws(end).cost,opts.rho,runtime);
end

% Update penalty parameter
if opts.adaptive; [ws, opts] = updatePenaltyParameter(iter, ws, opts); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ws, opts] = updatePenaltyParameter(iter, ws, opts)
% Update penalty parameter and take into account rescaling of the dual variables
persistent itPinf itDinf
if iter == 1 || isempty(itPinf)
    itPinf = 0; % number of iterations for which pinf/dinf <= eta
    itDinf = 0; % number of iterations for which pinf/dinf > eta
end
update = false;
resRat = ws(end).pres/ws(end).dres;
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
    parfor i = 1:opts.noCliques
        ws(i).eta = scaling .* ws(i).eta;
        ws(i).xi = scaling .* ws(i).xi;
    end
    opts.rho = newRho;
    opts.lambdaComp = 1-opts.lambda;
    opts.lambdaOverRho = opts.lambda / opts.rho;
    opts.lambdaCompOverRho = opts.lambdaComp / opts.rho;
end
end