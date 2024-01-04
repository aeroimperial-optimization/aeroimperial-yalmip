function [header,myline1,myline2] = printHeader()

% Print headers for solver
myline1 = [repmat('=',1,64),'\n'];
myline2 = [repmat('-',1,64),'\n'];
header  = ' iter |   pres   |   dres   |    cost    |   rho    | time (s) |\n';
