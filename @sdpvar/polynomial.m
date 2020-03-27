function [p,c,v] = polynomial(x,dmax,dmin,symmetries)
%POLYNOMIAL Creates parameterized polynomial
%
% [p,c,v] = polynomial(x,dmax,dmin)
%
% POLYNOMIAL is a quick way to define a parameterized polynomial p=c'*v,
% with all monomials of dmin <= degree(p,x) <= dmax. The coefficients in
% the polynomial are c while v is the monomial basis.
%
% Example:
%
% Paramterized quartic
%  x = sdpvar(2,1);
%  p = polynomial(x,4);
%
% See also MONOLIST, COEFFICIENTS

if (length(dmax) > 1) && (length(dmax) ~= length(x))
    error('Dimension mismatch: The second argument should be the max degree for each variable, or a sclar');
end

if nargin > 2
    if (length(dmin) > 1) && (length(dmin) ~= length(x))
        error('Dimension mismatch: The second argument should be the max degree for each variable, or a sclar');
    end
end

if any(dmax < 0)
    error('Only non-negative polynomial degrees possible')
end

if nargin<3
    dmin = 0;
end

if nargin < 4
    symmetries = [];
end

if isempty(dmin)
    dmin = 0;
end

if any(dmin > dmax)
    error('Third argument (dmin) should not be larger than second argument (dmax)');
end

if any(dmin < 0)
    error('Only non-negative polynomial degrees possible')
end

% multipartite?
if length(dmax)>1 && any(dmax(1)~=dmax)
    % symmetries will be ignored
    v = monolist(x,dmax,dmin,symmetries);
    c = sdpvar(length(v),1);
    p = c'*v;
    return
else
    % slightly faster?
    dmax =  dmax(1);
    powers = monpowers(length(x),dmax,symmetries);
    powers = powers(sum(powers,2)>=dmin,:);
    v = recovermonoms(powers,x);
    if dmin <= dmax && dmin>0
        s = nchoosek(length(x) + dmin-1,dmin-1);
        v = extsubsref(v,s+1:length(v));
    end
    c = sdpvar(length(v),1);
    p = c'*v;
end