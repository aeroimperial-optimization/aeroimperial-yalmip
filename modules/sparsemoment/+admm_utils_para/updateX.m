function ws = updateX(ws, opts, iter)

% UPDATEX
% Update the variables s{i}. For each index i, must solve a linear system
% We cache the factorization of the matrix at the first iteration

tstart = tic;
y = ws(end).y;
lambdaCompOverRho = opts.lambdaCompOverRho;
parfor i = 1:opts.noCliques
    % Initalize if needed
    if iter==1
        ws(i) = KKTfactor(ws(i));
    end
    
    % Right-hand side vector
    rhs = lambdaCompOverRho .* ws(i).bs;
    rhs = rhs + ws(i).At.' * (ws(i).c - ws(i).z + ws(i).eta);
    rhs = rhs + ws(i).E.' * (ws(i).d - ws(i).F*y + ws(i).xi);
    
    % Solve
    ws(i).s = KKTsolve(ws(i).s, rhs, ws(i).factor);
end
% Finished
ws(end).time.updateX = ws(end).time.updateX + toc(tstart);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ws = KKTfactor(ws)
H = ws.E.' * ws.E + ws.At.' * ws.At;
[ws.factor.L, ~, ws.factor.p] = chol(H,'lower','vector');
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