function all_monoms = monpowers(n,d,symmetries)

if nargin < 3; symmetries = []; end


% =========================================================================== %
% Original YALMIP code
% --------------------------------------------------------------------------- %
% if max(d)==0
%     all_monoms = [];
% else
%     all_monoms  = fliplr(eye(n));
%     last_monoms = all_monoms;
% end
% for degrees = 1:1:d-1
%     new_last_monoms = [];
%     for variable = 1:n
%         temp = last_monoms;
%         temp(:,variable) = temp(:,variable)+1;
%         new_last_monoms = [new_last_monoms;temp];
%      %   all_monoms = [all_monoms;temp];
%     end
%     last_monoms = unique(new_last_monoms,'rows');
%     all_monoms = [all_monoms;last_monoms];
%     %all_monoms = unique(all_monoms,'rows');
% end
% all_monoms = [zeros(1,n);all_monoms];
% all_monoms = fliplr(all_monoms);
% =========================================================================== %



% =========================================================================== %
% Same as above but preallocate for speed
% By Giovanni Fantuzzi
% --------------------------------------------------------------------------- %
if max(d)==0
    all_monoms = zeros(1,n);
    
else
    
    % Indexing
    global_ind = zeros(d+1,1);
    for i = 2:d+1
        global_ind(i) = global_ind(i-1) + nchoosek(n-2+i,i-1);
    end
    global_ind = global_ind + 1;
    
    % Initialize
    all_monoms = zeros(nchoosek(n+d,d),n);
    last_monoms  = fliplr(eye(n));
    all_monoms(2:n+1,:) = last_monoms;
    ptr = 1:n;
    
    % Loop
    for degrees = 1:1:d-1
        counter = global_ind(degrees+1);
        for i = 1:n
            new_counter = counter+ptr(i);
            all_monoms(counter+1:new_counter,:) = last_monoms(1:ptr(i),:);
            all_monoms(counter+1:new_counter,n-i+1) = last_monoms(1:ptr(i),n-i+1)+1;
            counter = new_counter;
        end
        ptr = [1, ptr(2:n)+cumsum(ptr(1:n-1))];
        last_monoms = all_monoms(global_ind(degrees+1)+1:global_ind(degrees+2),:);
    end
    
    % Flip to conclude (keep same order as original in YALMIP)
    all_monoms = fliplr(all_monoms);
end
% =========================================================================== %


% Finally, apply requested sign symmetries
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