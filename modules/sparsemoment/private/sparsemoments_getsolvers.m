function [solver, problem] = sparsemoments_getsolvers(options)


%% persistent variables (speed)
persistent CACHED_SOLVERS
persistent EXISTTIME
persistent NCHECKS
persistent allsolvers

% *************************************************************************
%% LOOK FOR AVAILABLE SOLVERS
% Finding solvers can be very slow on some systems. To alleviate this
% problem, YALMIP can cache the list of available solvers.
% *************************************************************************
if (options.cachesolvers==0) || isempty(CACHED_SOLVERS)
    getsolvertime = clock;
    [solvers,~,allsolvers] = getavailablesolvers(0,options);
    getsolvertime = etime(clock,getsolvertime);
    % CODE TO INFORM USERS ABOUT SLOW NETWORKS!
    if isempty(EXISTTIME)
        EXISTTIME = getsolvertime;
        NCHECKS = 1;
    else
        EXISTTIME = [EXISTTIME getsolvertime];
        NCHECKS = NCHECKS + 1;
    end
    if (options.cachesolvers==0)
        if ((NCHECKS >= 3 && (sum(EXISTTIME)/NCHECKS > 1)) || EXISTTIME(end)>2)
            if warningon
                info = 'Warning: YALMIP has detected that your drive or network is unusually slow.\nThis causes a severe delay in OPTIMIZE when I try to find available solvers.\nTo avoid this, use the options CACHESOLVERS in SDPSETTINGS.\nSee the FAQ for more information.\n';
                fprintf(info);
            end
        end
    end
    if length(EXISTTIME) > 5
        EXISTTIME = EXISTTIME(end-4:end);
        NCHECKS = 5;
    end
    CACHED_SOLVERS = solvers;
else
    solvers = CACHED_SOLVERS;
end



% *************************************************************************
%% NO SOLVER AVAILABLE
% *************************************************************************
if isempty(solvers)
    error('No solver available?')
end



% ***********************************************
%% SET UP PROBLEM CATEGORY
% A standard conic problem with lienar equalities, inequalities, SDP and
% SOCP constraints (don't actually have SOCP?)
% ***********************************************
problem.objective.linear = 1;
problem.objective.quadratic.convex = 0;
problem.objective.quadratic.nonconvex = 0;
problem.objective.quadratic.nonnegative = 0;
problem.objective.polynomial = 0;
problem.objective.maxdet.convex = 0;
problem.objective.maxdet.nonconvex = 0;
problem.objective.sigmonial = 0;

problem.constraint.equalities.linear     = 1;
problem.constraint.equalities.quadratic  = 0;
problem.constraint.equalities.polynomial = 0;
problem.constraint.equalities.sigmonial  = 0;
problem.constraint.equalities.multiterm = 0;

problem.constraint.inequalities.elementwise.linear = 1;
problem.constraint.inequalities.elementwise.quadratic.convex = 0;
problem.constraint.inequalities.elementwise.quadratic.nonconvex = 0;
problem.constraint.inequalities.elementwise.sigmonial = 0;
problem.constraint.inequalities.elementwise.polynomial = 0;

problem.constraint.inequalities.semidefinite.linear = 1;
problem.constraint.inequalities.semidefinite.quadratic = 0;
problem.constraint.inequalities.semidefinite.polynomial = 0;
problem.constraint.inequalities.semidefinite.sigmonial = 0;

problem.constraint.inequalities.rank = 0;

problem.constraint.inequalities.secondordercone.linear = 0;
problem.constraint.inequalities.secondordercone.nonlinear = 0;
problem.constraint.inequalities.rotatedsecondordercone.linear = 0;
problem.constraint.inequalities.rotatedsecondordercone.nonlinear = 0;
problem.constraint.inequalities.powercone = 0;

problem.constraint.complementarity.variable = 0;
problem.constraint.complementarity.linear = 0;
problem.constraint.complementarity.nonlinear = 0;

problem.constraint.integer = 0;
problem.constraint.binary = 0;
problem.constraint.semicont = 0;
problem.constraint.sos1 = 0;
problem.constraint.sos2 = 0;

problem.complex = 0;
problem.parametric = 0;
problem.interval = 0;
problem.evaluation = 0;
problem.exponentialcone = 0;
problem.gppossible = 0;


% ***********************************************
%% GET SOLVER
% Call Yalmip's native function
% ***********************************************
[solver,problem] = selectsolver(options,problem,solvers,0,allsolvers);
end


% Useful stuff
function yesno = warningon
s = warning;
yesno = isequal(s(1).state,'on');
end