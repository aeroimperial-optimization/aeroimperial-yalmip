function [newton_m2,N_unique,newton_m2_unique] = monomialproducts(N,n)
%MONOMIALPRODUCTS  Internal function used for monomial reduction

% Author Johan L�fberg
% $Id: monomialproducts.m,v 1.1 2006-03-30 13:56:54 joloef Exp $

% Exponents in squared monomials

N_unique = [];
for i = 1:size(N,1)
    [nN,mN] = size(N{i});
    newton_m2{i} = zeros(nN^2,mN+2);
    shift = 0;
    for j = 1:nN
        newton_m2{i}(shift+1:shift+nN,:) = [(1:nN)' j.*ones(nN,1) N{i}(1:nN,:)+N{i}(j,:)];
        shift = shift+nN;
    end    
    % Whoops, double copies of diagonal (we want double copies of non-diagonals though)
    if isempty(newton_m2{i})
        newton_m2_unique{i} = [];
    else
        [dummy,j,dummy2] = uniquesafe(newton_m2{i}(:,1:2),'rows');
        newton_m2{i} = newton_m2{i}(j,:);
        % Extract unique monomial products
        [dummy,j,dummy2] = uniquesafe(newton_m2{i}(:,3:end),'rows');
        newton_m2_unique{i} = newton_m2{i}(j,:);
    end
    N_unique = [N_unique;newton_m2_unique{i}];   
end
if ~isempty(N_unique)
    [dummy,j,dummy2] = uniquesafe(N_unique(:,3:end),'rows');
    N_unique = N_unique(j,:);
end
