function all_monoms = monpowers(n,d,symmetries)

if nargin < 3; symmetries = []; end

if max(d)==0
    all_monoms = [];
else
    all_monoms  = fliplr(eye(n));
    last_monoms = all_monoms;
end

for degrees = 1:1:d-1
    new_last_monoms = [];
    for variable = 1:n
        temp = last_monoms;
        temp(:,variable) = temp(:,variable)+1;
        new_last_monoms = [new_last_monoms;temp];
        %   all_monoms = [all_monoms;temp];
    end
    last_monoms = unique(new_last_monoms,'rows');
    all_monoms = [all_monoms;last_monoms];
    %all_monoms = unique(all_monoms,'rows');
end

all_monoms = [zeros(1,n);all_monoms];
all_monoms = fliplr(all_monoms);

if ~isempty(symmetries)
    all_monoms = applySymmetries(all_monoms,symmetries);
end

end

% ----------------------------------------------------------------------- %
function powers = applySymmetries(powers,symmetries)
% find rows of powers that represent monomials that are invariant under
% all symmetries specified by the user
H = rem(powers*symmetries,2)==0;
powers = powers(all(H,2),:);
end