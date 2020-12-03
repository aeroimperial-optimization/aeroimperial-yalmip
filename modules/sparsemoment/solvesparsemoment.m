function [pstar, y, all_moments, sol, model] = solvesparsemoment(x, p, h, g, omega, mass, options, cliques)

% [pstar,y,exponents,sol,model] = solvesparsemoment(x,p,h,g,omega,mass,options,cliques)
%       solves the moment relaxation of the standard-form polynomial 
%       optimization problem (POP)
%
%       min_x p(x) s.t. h_i(x)=0, g_j(x)>=0
%
%       by exploiting correlative sparsity. 
%
%       NOTES: 
%       (1) This is research code. It is not optimized, not guaranteed to
%           work, and probably buggy. Use at your own risk!
%       (2) This code does NOT exploit symmetry
%
%       Required inputs
%       ----------------
%       * x: the optimization variables of the POP (sdpvar vector)
%       * p: the polynomial objective to be MINIMIZED
%       * h: a vector of polynomial equality constraints, h_i(x)=0
%       * g: a vector of polynomial inequality constraints, g_i(x)>=0
%       * omega: the degree of the relaxation. The max degree of monomials
%                considered in the moment relaxation and associated SOS
%                program will be 2*omega.
%
%       Optional inputs
%       ---------------
%       * MassValue: total mass of the the moment-generating measure. The default
%                    value is 1, but setting a different positive value may
%                    improve the numerical behaviour (it amounts to problem
%                    scaling)
%       * options: yalmip options created with sdpsettings()
%       * cliques: sets of variables specifying the problem's sparsity
%                  structure. If not provided, the sparsity is detected
%                  automatically using a chordal extension of the
%                  correlative sparsity graph (see Waki et al, 2006)
%
%       Outputs
%       -------
%       * pstar: optimal value of the moment-SDP relaxation. it is a lower
%                bound on the optimal value of the POP
%       * y: the vector of moments solving the moment-SDP relaxation
%       * sol: the solver's output (if required)
%       * model: the conic problem model (in SeDuMi format)
%       * exponents: the exponents of the monomials modelled by y in the
%                    moment-SDP relaxation. For POPs with a unique solution
%                    and whose sparsity satisfies the running intersection
%                    property, one has
%               
%                    y(i) -> x(1)^exponents(i,1) * ... * x(n)^exponents(i,n)
%
%                    as the relaxation parameter omega tends to infinity.
%       
%
% Giovanni Fantuzzi
% 02 Dec 2019

% Inputs correct?
if nargin < 5; omega = []; end                  % No relaxation order? Set as empty
if nargin < 6; mass = 1; end                    % No mass specified? Use 1
if nargin < 7; options = sdpsettings; end       % No options? use yalmip default
if nargin < 8; cliques = []; end                % No cliques? set to empty
if isempty(options); options = sdpsettings; end % Empty options? use yalmip default
if isempty(mass); mass = 1; end                 % Empty mass? set to 1

% Compile moment problem
[At,b,c,K,isMomentMatrix,CD,PROJ,bshift,y0,b0,setuptime,all_moments,gramMonomials] = ...
    compilesparsemoment(x, p, h, g, omega, mass, options, cliques);

% Finally, solve and set outputs. Try to use generic yalmip format to 
% interface easily with user-preferred choice. Tested with:
% sedumi, mosek, sdpt3, cdcs, scs 
% Could be buggy with other solvers.
% NOTE: ensure we save solver outputs to get the solution
if options.verbose; disp('Solving...'); end
options.savesolveroutput = 1;
[solver,problemClass] = sparsemoments_getsolvers(options);
interfacedata = sparsemoments_interfacedata(At,b,c,K,options,solver,problemClass);
try
    eval(['output = ' solver.call '(interfacedata);']);
catch
    error('Woops, something went wrong in the solver!')
end

% Get solution, scale y by the zero-th moment and get optimal value of the 
% moments-SDP relaxation.
y = y0 + PROJ*output.Primal;
y = y./mass;
pstar = (b.'*output.Primal + bshift)/mass + b0;

% Output solver solution if desired
if nargout > 3
    sol = output.solveroutput; 
    sol.setupTime = setuptime;
    sol.reducedMoments = output.Primal;
    sol.momentMatrices = recoverMomentMatrices(output.Primal, At, c, K, isMomentMatrix);
    sol.gramMonomials = gramMonomials;
    sol.cliques = CD.cliques;
end

% Output model if needed
% Change sign of b here to have dual-standard-form problem:
% max b'*y    s.t.    c - At*y \in K
if nargout > 4
    model.At = At;
    model.b = -b;
    model.c = c;
    model.K = K;
end