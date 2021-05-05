function m = monolistcoeff_fast(n,d,summax)
%MONOLISTCOEFF_FAST  Internal function used in SOS programs
%Make exponents of monomials in fast and memory efficient way
%NOTE: the exponents are ordered in a random way

z = sparse((0:d).');
E = ones(length(z),1);
m = z;
for i = 2:n
    m = [kron(m,E) kron(ones(size(m,1),1),z)];
    m = unique(m,'rows');
    m = m(sum(m,2)<=summax,:);
end