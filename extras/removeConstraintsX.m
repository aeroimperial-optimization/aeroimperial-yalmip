function [A,b,x0,PROJ,done] = removeConstraintsX(A,b)

% Reduce the set of equations A*x=b by expressing the solution as
% x=x0+PROJ*y, where y solves A*y=b for modified A and b. Use constraints 
% that involve at most TWO variables preserve sparsity

% Initalize projection matrix for moments
num_x = size(A,2);
PROJ = speye(num_x);
x0 = sparse(num_x,1);
done = true;

% Find constraints that only depend on two moments
twovar = sum( spones(A), 2 )<=2;
while any(twovar)
    A_twomoments = A(twovar,:);
    b_twomoments = b(twovar);
    A_others = A(~twovar,:);
    b_others = b(~twovar);
    done = false; % we can eliminate some variables
    
    % Dependent constraints? Ignore them...
    indep_idx = lirows([b_twomoments, A_twomoments]);
    if ~all(indep_idx)
        b_twomoments = b_twomoments(indep_idx);
        A_twomoments = A_twomoments(indep_idx,:);
    end
    
    % Get a tolerance for cleaning
    [m,n] = size(A_twomoments);
    TOL = max(m,n)*eps(class(A_twomoments))*norm(A_twomoments,inf);

    % Now solve for x = y+N*z, where y is known and N satisfies 
    % A_twomoments*N = 0. We hope that N has only one entry on each row.
    y = A_twomoments\b_twomoments;
    N = fastnull((A_twomoments),'r');
    P = sum(spones(N),2);
    if ~all(P<=1)
        warning('Cannot eliminate variables in a trivial way!')
        return
    end
    b = b_others - A_others*y;
    A = A_others*N;
    A = cleanzeros(A, TOL);
    b = cleanzeros(b, TOL);
    
    % Eliminate zero rows?
    zeroRows = sum(spones(A),2)==0;
    if any(zeroRows)
        if any(b(zeroRows))
            error('Inconsistent constraints: requested C = 0 with C~=0!')
        else
            A = A(~zeroRows,:);
            b = b(~zeroRows);
        end
    end
    
    % Update projection
    x0 = x0 + PROJ*y;
    PROJ = PROJ*N;
    
    % Did we reveal anything else?
    twovar = sum( spones(A), 2)<=2;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning function for sparse matrices
function A = cleanzeros(A, TOL)
    % Clean values smaller than a tolerance
    if issparse(A)
        [m,n] = size(A);
        [iA,jA,vA] = find(A);
        idx = abs(vA)>=TOL;
        A = sparse(iA(idx),jA(idx),vA(idx),m,n);
        
    else
        A(abs(A)<TOL)=0;
    end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = lirows(X,tol)
% Find indices of a linearly independent set of rows of a given matrix X
%
%    idx=lirows(X)
%
% inputs:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
% outputs:
%
% idx:  Logical values of the independent rows of X

% X has no non-zeros
if ~nnz(X)
    idx=[];
    return
end

if nargin<2
    tol=1e-10;
end

% [Q,R,E] = qr(full(X.'),0);   % full for compatibility with old MATLAB
[Q,R,E] = qr(X.',0);
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = abs(R(1));
end

r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
idx = sort(E(1:r));
idx = ismember(1:size(X,1),idx);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = fastnull(A,how)
%NULL   Null space.
%   Z = NULL(A) is an orthonormal basis for the null space of A obtained
%   from the singular value decomposition.  That is,  A*Z has negligible
%   elements, size(Z,2) is the nullity of A, and Z'*Z = I.
%
%   Z = NULL(A,'r') is a "rational" basis for the null space obtained
%   from the reduced row echelon form.  A*Z is zero, size(Z,2) is an
%   estimate for the nullity of A, and, if A is a small matrix with
%   integer elements, the elements of R are ratios of small integers.
%
%   The orthonormal basis is preferable numerically, while the rational
%   basis may be preferable pedagogically.
%
%   Example:
%
%       A =
%
%           1     2     3
%           1     2     3
%           1     2     3
%
%       Z = null(A);
%
%       Computing the 1-norm of the matrix A*Z will be
%       within a small tolerance
%
%       norm(A*Z,1)< 1e-12
%       ans =
%
%          1
%
%       null(A,'r') =
%
%          -2    -3
%           1     0
%           0     1
%
%   Class support for input A:
%      float: double, single
%
%   See also SVD, ORTH, RANK, RREF.

%   Copyright 1984-2017 The MathWorks, Inc.

[m,n] = size(A);
if nargin > 1 && (isequal(how,'r') || isequal(how,"r"))
    
    % Rational basis
    [R,pivcol] = fastrref(A);
    r = length(pivcol);
    nopiv = 1:n;
    nopiv(pivcol) = [];
    if ~issparse(A)
        % code for dense matrices
        Z = zeros(n,n-r,class(A));
        if n > r
            Z(nopiv,:) = eye(n-r,n-r,class(A));
            if r > 0
                Z(pivcol,:) = -R(1:r,nopiv);
            end
        end
    else
        % code for sparse matrices
        iZ = []; jZ = []; vZ = [];
        if n > r
            % ZZ(nopiv,:) = speye(n-r,n-r);
            iZ = [iZ; nopiv(:)];
            jZ = [jZ; (1:n-r)'];
            vZ = [vZ; ones(n-r,1)];
            if r > 0
                % ZZ(pivcol,:) = -R(1:r,nopiv);
                iZ = [iZ; repmat(pivcol(:),n-r,1)];
                col_idx = repmat(1:n-r,length(pivcol),1);
                jZ = [jZ; col_idx(:)];
                vals = -R(1:r,nopiv);
                vZ = [vZ; vals(:)];
            end
        end
        Z = sparse(iZ,jZ,vZ,n,n-r);
    end
    
else
    
    % Orthonormal basis
    [~,S,V] = svd(A,0);
    if isempty(A)
        Z = V;
    else
        if m == 1
            s = S(1);
        else
            s = diag(S);
        end
        tol = max(m,n) * eps(max(s));
        r = sum(s > tol);
        Z = V(:,r+1:n);
    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, jb] = fastrref(A, tol)
%RREF   Reduced row echelon form.
%   R = RREF(A) produces the reduced row echelon form of A.
%
%   [R,jb] = RREF(A) also returns a vector, jb, so that:
%       r = length(jb) is this algorithm's idea of the rank of A,
%       x(jb) are the bound variables in a linear system, Ax = b,
%       A(:,jb) is a basis for the range of A,
%       R(1:r,jb) is the r-by-r identity matrix.
%
%   [R,jb] = RREF(A,TOL) uses the given tolerance in the rank tests.
%
%   Roundoff errors may cause this algorithm to compute a different
%   value for the rank than RANK, ORTH and NULL.
%
%   Class support for input A:
%      float: double, single
%
%   See also RANK, ORTH, NULL, QR, SVD.

%   Copyright 1984-2017 The MathWorks, Inc.

useqr = 0;
[m,n] = size(A);

% % Does it appear that elements of A are ratios of small integers?
% [num, den] = rat(A);
% rats = isequal(A, num./den);

% Compute the default tolerance if none was provided.
if (nargin < 2)
    tol = max(m,n)*eps(class(A))*norm(A,inf);
end

if ~useqr
% Loop over the entire matrix.
i = 1;
j = 1;
jb = zeros(1,0);
while i <= m && j <= n
    % Find value and index of largest element in the remainder of column j.
    [p, k] = max(abs(A(i:m,j)));
    k = k+i-1;
    if p <= tol
        % The column is negligible, zero it out.
        A(i:m,j) = 0;
        j = j + 1;
    else
        % Remember column index
        jb = [jb j]; %#ok<AGROW>
        % Swap i-th and k-th rows.
        A([i k],j:n) = A([k i],j:n);
        % Divide the pivot row by the pivot element.
        A(i,j:n) = A(i,j:n)./A(i,j);
        % Subtract multiples of the pivot row from all the other rows.
        % -------------------------------------------------
        % MATLAB's VERSION: SLOW LOOP
        %       for k = [1:i-1 i+1:m]
        %          A(k,j:n) = A(k,j:n) - A(k,j).*A(i,j:n);
        %       end
        % -------------------------------------------------
        % Vectorized loop using "outer product" capabilities
        k = [1:i-1 i+1:m];
        A(k,j:n) = A(k,j:n) - A(k,j).*A(i,j:n);
        i = i + 1;
        j = j + 1;
    end
end
end

% % Return "rational" numbers if appropriate.
% if rats
%     [num, den] = rat(A);
%     A = num./den;
% end
end