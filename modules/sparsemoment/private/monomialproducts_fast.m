function [BBt,N_unique] = monomialproducts_fast(B)
%MONOMIALPRODUCTS  Internal function used for monomial reduction
%Try to work with sparse arrays
%Work with a single input (not the same as YALMIP!)

% Make things fast
if ~issparse(B); B = sparse(B); end
[m,n] = size(B);
x = rand(n,1);
shift = 0;

% Operate
A = (1:m).';
I = []; J = []; V = [];
for j = 1:m
    temp = [A, repmat(j,m,1), B+B(j,:)];
    [ii,jj,vv] = find(temp);
    I = [I; ii+shift];
    J = [J; jj];
    V = [V; vv];
    shift = shift+m;
end
BBt = sparse(I, J, V, m^2, n+2);
    
% Unique squared moments
[~, ind] = uniquetol( full(BBt(:,3:end)*x) );
N_unique = BBt(ind,:);