function interfacedata = sparsemoments_interfacedata(At,b,c,K,options,solver,problemClass)

% *************************************************************************
%% GENERAL DATA EXCHANGE WITH SOLVER
% *************************************************************************
nvars = size(At,2);
interfacedata.F_struc = [c, -At];
interfacedata.c = -b;
interfacedata.Q = sparse(nvars,nvars);
interfacedata.f = 0;
interfacedata.K = K;
interfacedata.lb = -inf(nvars,1);
interfacedata.ub =  inf(nvars,1);
interfacedata.x0 = [];
interfacedata.options = options;
interfacedata.solver  = solver;
interfacedata.monomtable = [];
interfacedata.variabletype = [];
interfacedata.integer_variables   = [];
interfacedata.binary_variables    = [];
interfacedata.semicont_variables    = [];
interfacedata.semibounds = [];
interfacedata.uncertain_variables = [];
interfacedata.parametric_variables= [];
interfacedata.extended_variables  = [];
interfacedata.aux_variables  = [];
interfacedata.used_variables      = [];
interfacedata.lowrankdetails = [];
interfacedata.problemclass = [];
interfacedata.KCut = [];
interfacedata.getsolvertime = 1;
% Data to be able to recover duals when model is reduced
interfacedata.oldF_struc = [];
interfacedata.oldc = [];
interfacedata.oldK = [];
interfacedata.factors = [];
interfacedata.Fremoved = [];
interfacedata.evalMap = [];
interfacedata.evalVariables = [];
interfacedata.evaluation_scheme = [];
interfacedata.lift = [];
if strcmpi(solver.tag,'bmibnb')
    interfacedata.equalitypresolved = 0;
    interfacedata.presolveequalities = 1;
else
    interfacedata.equalitypresolved = 1;
    interfacedata.presolveequalities = 1;
end
interfacedata.ProblemClass = problemClass;
interfacedata.dualized = 0;

end
