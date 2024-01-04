function sol = updateX(prog, sol, opts, iter)

% UPDATEX
% Update the variables s{i}. For each index i, must solve a linear system
% We cache the factorization of the matrix at the first iteration

tstart = tic;
persistent factors
for i = 1:opts.noCliques
    % Initalize if needed
    if iter==1 || isempty(factors)
        tmp = KKTfactor(i, prog, opts);
        factors(i).L = tmp.L;
        factors(i).p = tmp.p;
    end
    
    % Right-hand side vector
    rhs = opts.lambdaCompOverRho .* prog.b.s{i};
    rhs = rhs + prog.At{i}.' * (prog.c{i} - sol.z{i} + sol.eta{i});
    rhs = rhs + prog.E{i}.' * (prog.d{i} - prog.F{i}*sol.y + sol.xi{i});
    
    % Solve
    sol.s{i} = KKTsolve(sol.s{i}, rhs, factors(i));
end
% Finished
sol.time.updateX = sol.time.updateX + toc(tstart);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function factor = KKTfactor(i, prog, opts)
H = prog.E{i}.' * prog.E{i} + prog.At{i}.' * prog.At{i};
[factor.L, ~, factor.p] = chol(H,'lower','vector');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = KKTsolve(s, rhs, f)
% Solve the KKT system for clique i
% TODO: Implement other method if Cholesky fails or is unsuitable!

% Use built-in solvers?
persistent useBuiltin
if(isempty(useBuiltin))
    %default to look for CSparse code
    useBuiltin = ~exist(['.',filesep,'cs_ltsolve.' mexext],'file');
    useBuiltin = useBuiltin | ~exist(['.',filesep,'cs_lsolve.' mexext],'file');
end

% The actual operation
if(useBuiltin)
    %Native matlab version (slow)
    s(f.p) = f.L'\(f.L\rhs(f.p,:));
else
    %Csparse version (avoids transpose)
    s(f.p) = cs_ltsolve(f.L,cs_lsolve(f.L,rhs(f.p,:)));
end

end