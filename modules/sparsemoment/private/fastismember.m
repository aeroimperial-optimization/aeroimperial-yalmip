function AsubsetB = fastismember(A,B)

% Fast way to check if A is a subset of B, where A and B are sets of
% positive integers. Assume nonempty

% Intersect A and B
%if ~isempty(A)&&~isempty(B)
P = false(1, max(max(A),max(B)) ) ;
P(B) = true;
C = A(P(A));
nC = length(C);
nA = length(A);
AsubsetB = false;
if nA==nC
    AsubsetB = sum(A-C)==0;
end
%else
%    C = [];
%end