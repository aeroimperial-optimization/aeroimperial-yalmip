% sparse_moment_relaxation.m
% --------------------------
% A tutorial in the use of the sparse moment relaxation module. We solve a
% correlatively sparse polynomial optimization problem,
%
% min   x(1)^2 + x(2)^2 + x(3)^2 - x(1) - x(2) - x(3)
%
% s.t.  x(1) - x(2)^2 = 0
%       x(2)*x(3) = 0
%       1 - x(1)^2 - x(2)^2 >= 0
%       1 - x(2)^2 - x(3)^2 >= 0
%
% We call "solvesparsemoment" to set up and solve the sparse moment-SOS and
% attempt to extract the optimal solution by reading the degree-1 moments.
% NOTE: This is not guaranteed to work in general!

% Clean up
clear
yalmip clear

% Parameters
% omega: the relaxation parameter (use moments up to degree 2*omega)
omega = 2;

% Ser up the polynomial optimization problem
% p = objective
% h = vector of equality constraints
% g = vector of inequality constraints
x = sdpvar(3,1);
p = dot(x,x) - sum(x);
h = [x(1)-x(2)^2; x(2)*x(3)];
g = [1-x(1)^2-x(2)^2; 1-x(2)^2-x(3)^2];

% Call the sparse moment solver and attempt to read the optimal x from the
% degree-1 moments using the "extractmomentsolution" function (this should work
% for large enough relaxation order omega if the optimal x is unique)
[pstar, y, exponents, sol, mod] = solvesparsemoment(x,p,h,g,omega);
x = extractlinearmoments(y, exponents);