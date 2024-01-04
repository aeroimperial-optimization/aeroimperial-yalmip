function [sol, opts] = makeVariables(prog, opts, sol)

% MAKEVARIABLES
% Initialize the variables for the solver
% TODO: allow for user-specified initial variables

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
        [sol, opts] = makeVariables(prog, opts, []);
    end
else
    % Initialize
    sol.y = zeros(opts.m.y,1);
    sol.s = cell(opts.noCliques,1);
    sol.z = cell(opts.noCliques,1);
    sol.eta = cell(opts.noCliques,1);
    sol.xi = cell(opts.noCliques,1);
    sol.cost = -opts.lambda * ( prog.b.y.' * sol.y );
    % Residuals (zero, will update pres below)
    sol.pres = 0;
    sol.dres = 0;
    % Clique variables: free (s{i}) and conic (z{i})
    % Dual variables: eta{i} and xi{i}
    % Everything is initalized to zero (on the boundary of the cone!)
    for i = 1:opts.noCliques
        zSize = prog.K(i).f + prog.K(i).l + prog.K(i).q + sum(prog.K(i).s.^2);
        sol.z{i} = zeros(zSize,1);
        sol.s{i} = zeros(opts.m.s(i),1);
        sol.eta{i} = zeros(opts.n.s(i),1);
        sol.xi{i} = zeros(length(prog.d{i}),1);
        sol.cost = sol.cost - opts.lambdaComp * ( prog.b.s{i}.' * sol.s{i} );
        sol.pres = max(sol.pres, norm(prog.c{i}-prog.At{i}*sol.s{i}-sol.z{i},'Inf'));
        sol.pres = max(sol.pres, norm(prog.d{i}-prog.E{i}*sol.s{i}-prog.F{i}*sol.y,'Inf'));
    end
end

% Times
sol.time.setupTime = 0;
sol.time.ADMMtime = 0;
sol.time.totalTime = 0;
sol.time.updateX  = 0;
sol.time.updateY  = 0;
sol.time.updateZ  = 0;

% Final penalty parameter
sol.penalty = opts.rho;

