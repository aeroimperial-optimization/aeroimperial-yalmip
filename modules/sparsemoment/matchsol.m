function x = matchsol(cliques,xsol)

% x = matchsol(cliques,xsol)
%
% Try to match up extracted solution across cliques. Stupid approach, but
% straightforward.

% Parameters
tol = 1e-2;

% Size of xsol not matching?
sz = cellfun(@(X) size(X,2), xsol);
if any(sz~=sz(1))
    error('Different number of minimizers in different cliques! Aborting.')
end

% Loop
for i = 1:length(cliques)
    for j = 1:length(cliques)
        [common, IA, IB] = intersect(cliques{i},cliques{j});
        if (i~=j) && ~isempty(common)
            % Overlapping cliques! check for matches
            for k = 1:size(xsol{i},2)
                err = abs( (xsol{i}(IA,k) - xsol{j}(IB,:))./(1+xsol{i}(IA)) )<tol;
                match = (sum(err)==length(common));
                if ~match
                    error('Inconsistent minimizers in different cliques! Aborting.')
                elseif match(k)~=1
                    % reorder
                    l = find(match,1,'first');
                    xsol{j}(:,[k l]) = xsol{j}(:,[l k]);
                end
            end
        end
    end
end

% Take average to account for small mismatch
cID = [cliques{:}].';
t = vertcat(xsol{:});
D = accumarray(cID,1);
x = zeros(length(D),sz(1));
for i = 1:sz(1)
   x(:,i) = accumarray(cID,t(:,i))./D;
end