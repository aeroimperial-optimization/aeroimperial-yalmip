function x = extractlinearmoments(y, exponents)

% x = extractlinearmoments(y, moments) reads the degree-1 monomials from a
% moment-SOS solution y. The values in y correspond to moments with
% exponents in the "exponents" argument.

islinear = (sum(exponents,2)==1);
[I,J] = find(exponents(islinear,:));
[~, order] = sort(I);
[~, order] = sort(J(order));
x = y(islinear);
x = x(order);