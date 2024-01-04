function [F_struc, K] = removeDependentEqualities(F_struc,K)
%removeDependentEqualities Internal function remove linearly dependent equality
% constraints

% Author Giovanni Fantuzzi, 1 Nov 2021

% Extract the equalities
A_equ = F_struc(1:K.f,2:end);
b_equ = -F_struc(1:K.f,1);

% Remove linear dependencies
remove = ~lirows([A_equ b_equ]);
F_struc(remove,:) = [];
K.f = sum(~remove);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proper identification of linearly dependent constraints
% GF on 1 Nov 2021
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