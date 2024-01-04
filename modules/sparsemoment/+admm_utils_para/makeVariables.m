function [ws, opts] = makeVariables(ws, opts, sol)

% MAKEVARIABLES
% Initialize the variables for the solver
% TODO: allow for user-specified initial variables

if nargin < 3; sol = []; end
if ~isempty(sol)
    try
        % TODO: Check that sol structure is correct.
        %       For now, assume that it is...
        % Restart using last penalty parameter and update options to
        % warm-start correctly
        opts.rho = sol.penalty;
        opts.lambdaComp = 1-opts.lambda;
        opts.lambdaOverRho = opts.lambda / opts.rho;
        opts.lambdaCompOverRho = opts.lambdaComp / opts.rho;
    catch
        warning('Woops, something went wrong with the variable initialization...')
        [ws, opts] = makeVariables(ws, opts);
    end
else
    % Initialize
    ws(end).y = zeros(opts.m.y,1);
    ws(end).cost = -opts.lambda * ( ws(end).b.' * ws(end).y );
    [ws(1:opts.noCliques).s] = deal([]);
    [ws(1:opts.noCliques).z] = deal([]);
    [ws(1:opts.noCliques).eta] = deal([]);
    [ws(1:opts.noCliques).xi] = deal([]);
    [ws(1:opts.noCliques).cost] = deal([]);
    [ws(1:opts.noCliques).pres] = deal([]);
    [ws(1:opts.noCliques).dres] = deal([]);
    [ws(1:opts.noCliques).dresvec] = deal([]);
    [ws(1:opts.noCliques).factor] = deal([]);
    % Residuals (zero, will update pres below)
    [ws(1:opts.noCliques).pres] = deal(0);
    [ws(1:opts.noCliques).dres] = deal(0);
    % Clique variables: free (s) and conic (z)
    % Dual variables: eta and xi
    % Everything is initalized to zero (on the boundary of the cone!)
    y = ws(end).y;
    ms = opts.m.s;
    ns = opts.n.s;
    lambdaComp = opts.lambdaComp;
    for i = 1:opts.noCliques
        zSize = ws(i).K.f + ws(i).K.l + ws(i).K.q + sum(ws(i).K.s.^2);
        ws(i).z = zeros(zSize,1);
        ws(i).s = zeros(ms(i),1);
        ws(i).eta = zeros(ns(i),1);
        ws(i).xi = zeros(length(ws(i).d),1);
        ws(i).cost = lambdaComp * ( ws(i).bs.' * ws(i).s );
        ws(i).pres = norm(ws(i).c-ws(i).At*ws(i).s-ws(i).z,'Inf');
        ws(i).pres = max(ws(i).pres, norm(ws(i).d-ws(i).E*ws(i).s-ws(i).F*y,'Inf'));
    end
end

% Overall residuals
ws(end).pres = max([ws(1:opts.noCliques).pres]);
ws(end).dres = max([ws(1:opts.noCliques).dres]);

% Times
ws(end).time.setupTime = 0;
ws(end).time.ADMMtime = 0;
ws(end).time.totalTime = 0;
ws(end).time.updateX  = 0;
ws(end).time.updateY  = 0;
ws(end).time.updateZ  = 0;

% Final penalty parameter
ws(end).penalty = opts.rho;

