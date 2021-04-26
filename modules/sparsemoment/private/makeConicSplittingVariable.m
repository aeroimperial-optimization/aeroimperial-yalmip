function prog = makeConicSplittingVariable(At, b, c, K, isMomentMatrix, inclique, nMoments, options)

% Build the conic program with splitting of local to global variables:
%
%   min   -(bshift.y + b'*y)
%   s.t.   c{i} - At{i}*s{i} \in K(i)
%          d{i} - E{i}*s{i} - F{i}*y = 0
%
% The matrix F{i} should be such that \sum_i F{i}'*F{i} is diagonal and invertible
% The matrix E{i} should be such that \sum_i E{i}'*E{i} is diagonal and invertible?
% The cost is split so the objective can be replaced with
%
% -( k*bshift.y + (1-k)*bshift.s ) - (k*b'*y + (1-k)*\sum_i b.s{i}'*s{i})

% Initialize program
prog.At = [];
prog.b.y = b;
prog.b.s = [];
prog.c = [];
prog.K = K;
prog.d = [];
prog.E = [];
prog.F = [];
prog.bshift.y = 0;
prog.bshift.s = 0;
prog.yPROJ = speye(nMoments);
prog.y0 = sparse(nMoments,1);
prog.isMomentMatrix = isMomentMatrix;

% Total projector before elimination of y variables
% Needed to split the cost vector across cliques in the correct way
P = accumarray(vertcat(inclique{:}), 1, [nMoments, 1]);
bscaled = b./P;

% Eliminate moments from y
if options.sparsemoment.eliminateMoments
    prog_full = makeConicUniqueMoments(At, b, c, K, isMomentMatrix, options);
    prog.yPROJ = prog_full.PROJ;
    prog.y0 = prog_full.y0;
    prog.b.y = prog_full.b;
    prog.bshift.y = prog_full.bshift;
end

% Loop over cliques
nCliques = length(At);
for i = 1:nCliques
    % Do we have any 1-by-1 PSD cones? Move to linear cones
    % We assume moment matrices are not 1-by-1...
    islin = (prog.K(i).s==1);
    if any(islin)
        prog.K(i).l = prog.K(i).l + sum(islin);
        prog.K(i).s = prog.K(i).s(~islin);
        prog.isMomentMatrix{i} = prog.isMomentMatrix{i}(~islin);
        At(i).l = [At(i).l; At(i).s(islin)];
        c(i).l = [c(i).l; c(i).s(islin)];
        At(i).s = At(i).s(~islin);
        c(i).s = c(i).s(~islin);
    end
    % Vectorize
    At(i).f = vertcat(At(i).f{:}); c(i).f = vertcat(c(i).f{:});
    At(i).l = vertcat(At(i).l{:}); c(i).l = vertcat(c(i).l{:});
    At(i).s = vertcat(At(i).s{:}); c(i).s = vertcat(c(i).s{:});
    prog.At{i} = [At(i).f; At(i).l; At(i).s];
    prog.At{i} = prog.At{i}(:,inclique{i});
    prog.c{i} = [c(i).f; c(i).l; c(i).s];
    prog.b.s{i} = bscaled(inclique{i});
    % Projector to local moment variables
    sizeClique = length(inclique{i});
    if options.sparsemoment.eliminateMoments
        % Find constraints implied by elimination of global moments:
        % d{i} - E{i}*s{i} - F{i}*y = 0
        prog.F{i} = -prog_full.PROJ(inclique{i},:);
        keep = sum(spones(prog.F{i}), 2)>0;
        nrows = sum(keep);
        prog.F{i} = prog.F{i}(keep,:);
        prog.d{i} = prog_full.y0(inclique{i}(keep));
        prog.E{i} = sparse(1:nrows, find(keep), 1, nrows, sizeClique);
        % Add extra equalities and eliminate variables from s{i}
        numCnstr = sum(~keep);
        AA = [sparse(1:numCnstr, find(~keep), 1, numCnstr, sizeClique); prog.At{i}];
        cc = [prog_full.y0(inclique{i}(~keep)); prog.c{i}];
        prog.K(i).f = prog.K(i).f + numCnstr;
        [prog.At{i},prog.b.s{i},prog.c{i},prog.K(i),M,bshift,prog.s0{i}] = eliminateMatchingMoments(AA,prog.b.s{i},cc,prog.K(i));
        % prog.At, prog.b, prog.c and prog.K already updated inside the
        % function. Here, need to update prog.d, prog.E and prog.bshift
        prog.bshift.s = prog.bshift.s + bshift;
        prog.d{i} = prog.d{i} - prog.E{i} * prog.s0{i};
        prog.E{i} = prog.E{i} * M;
        % Keep unique constraints using a probabilistic approach (works
        % almost surely)
        [~,keep] = lirows([prog.d{i}, -prog.E{i}, -prog.F{i}]);
        prog.d{i} = prog.d{i}(keep);
        prog.E{i} = prog.E{i}(keep,:);
        prog.F{i} = prog.F{i}(keep,:);
    else
        prog.d{i} = sparse(sizeClique,1);
        prog.E{i} = speye(sizeClique);
        prog.F{i} = -sparse(1:sizeClique, inclique{i}, 1, sizeClique, nMoments);
    end
end

% End function
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Find linearly independent rows of matrix
function [Xsub,idx] = lirows(X,tol)
%Extract a linearly independent set of rows of a given matrix X
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted rows of X
% idx:  The indices (into X) of the extracted rows
X = X.';
if ~nnz(X) %X has no non-zeros and hence no independent columns
    Xsub=[]; idx=[];
    return
end
if nargin<2, tol=1e-10; end
[~, R, E] = qr(X,0);
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = R(1);
end
%Rank estimation
r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
idx=sort(E(1:r));
Xsub=X(:,idx);
Xsub = Xsub.';
end
