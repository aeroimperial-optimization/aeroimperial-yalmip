function [p,c,v] = multidegpoly(x,dmax,dmin,symmetries)
% MULTIDEGPOLY Creates parameterized polynomial with different degrees in
% different variables, such that max power of x{i} is dmax(i)
% NOTE: x here is either a cell array of variables, or a vector of variables
%
% [p,c,v] = multdegpoly(x,dmax,dmin,symmetries)
%
% MULTIDEGPOLY is a quick way to define a parameterized polynomial p=c'*v,
% whose monomials have powers of x{i} between dmin(i) and dmax(i). The 
% coefficients in the polynomial are c while v is the monomial basis.
% The function also allows the creation of sign-symmetric polynomials. The
% sign symmetries are specified by a 4th argument "symmetries", a length(x)-by-m
% matrix of 1s and 0s specifying m sign symmetries of the monomials in the list 
% (1s denote entries of x for which the sign is changed under the symmetry
% transformation).
%
% Example:
%
% Paramterized quartic
%  x = sdpvar(2,1);
%  p = polynomial({x(1), x(2)},[1 2],[1 1]);
%  sdisplay(p)
%       x(1)*c(1)+x(2)*c(2)+x(1)*x(2)*c(3)+x(2)^2*c(4)
%
% See also MONOLIST, COEFFICIENTS
% Written by G. Fantuzzi on 21 Oct 2020

if (length(dmax) > 1) && (length(dmax) ~= length(x))
    error('Dimension mismatch: The second argument should be the max degree for each variable, or a scalar');
end

if nargin > 2
    if (length(dmin) > 1) && (length(dmin) ~= length(x))
        error('Dimension mismatch: The second argument should be the max degree for each variable, or a scalar');
    end
    if length(dmin) ~= length(dmax)
        error('Dimension mismatch: Vectors dmax and dmin should have the same length');
    end
end

if nargin<3 || isempty(dmin)
    dmin = zeros(size(dmax));
end

if nargin < 4
    symmetries = [];
end

if any(dmax < 0) || any(dmin < 0)
    error('Only non-negative polynomial degrees possible')
end

if any(dmin > dmax)
    error('Third argument (dmin) should not be larger than second argument (dmax)');
end

% Convert to cell
if ~iscell(x)
    x = num2cell(x);
end

% Make monomial powers in a very ugly loop
exponents = monpowers(length(x{1}), dmax(1)); % max degree is dmax(i)
keep = sum(exponents, 2) >= dmin(1) ; % min degree is dmin(i)
exponents = exponents(keep, :);
xvar = x{1}(:);
for i = 2:length(x)
    v = monpowers(length(x{i}), dmax(i));
    keep = sum(v, 2) >= dmin(i) ;
    exponents = multipartiteMonomials(exponents, v(keep,:));
    xvar = [xvar; x{i}(:)];
end
% Apply symmetries
if ~isempty(symmetries)
    exponents = applySymmetries(exponents,symmetries);
end
% Construct polynomial
v = recovermonoms(exponents,xvar);
c = sdpvar(length(v),1);
p = c'*v;

%end main function
end


% ----------------------------------------------------------------------- %
function powers = applySymmetries(powers,symmetries)
% find rows of powers that represent monomials that are invariant under
% all symmetries specified by the user
H = rem(powers*symmetries,2)==0;
powers = powers(all(H,2),:);
end

% ----------------------------------------------------------------------- %
function Z = multipartiteMonomials(X,Y)
% make matrix of exponents for outer product of monomials in X and monomials in
% Y, which depend on DIFFERENT variables. Reorder using sortrows to obtain
% same order of monomials as created by yalmip's "monomialproducts". This
% is needed to build moment matrices more easily.
Z = kroncat(X,Y);
Z = sortrows(Z);
end

% ----------------------------------------------------------------------- %
% kronecker concatenation of matrices
function Z = kroncat(X,Y)
[rX,cX] = size(X);
[rY,cY] = size(Y);
Z = repmat(X,1,rY);
Z = reshape(Z',cX,rY*rX).';
Z = [Z, repmat(Y,rX,1)];
end
