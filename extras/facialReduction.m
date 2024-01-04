function [F_struc,b,K] = facialReduction(F_struc,b,K)

% Perform a stupid facial reduction for SDPs in standard primal form, trying to
% detect problems whose equality constraints force zero diagonal entries on PSD
% matrices. If found, we add constraints setting to zero the entire column
% of the corresponding PSD matrix, and repeat the procedure until
% convergence. Assume that the problem is specified in SeDuMi format.

% Project equalities (simple)
At = -F_struc(:,2:end);
c = F_struc(:,1);

[~,~,x0,PROJ,done] = removeConstraintsX(At.',b);

% Loop
while ~done
    
    % Construct the variable x "symbolically". It has a zero in position i
    % if and only if the equality constraints A*x=b force x(i)=0.
    x = x0 + spones(PROJ)*rand(size(PROJ,2),1);
    keep = true(size(At,1),1);
    shift = 0;
    shiftPSD = computePSDShift(K);
    
    % Remove zeros from free cone
    if isfield(K,'f') && ~isempty(K.f) && K.f>0
        idx = shift + (1:K.f);
        shift = shift + K.f;
        keep(idx) = x(idx)~=0;
        K.f = sum(keep(idx));
    end
    
    % Remove zeros from positive orthant
    if isfield(K,'l') && ~isempty(K.l)  && K.l>0
        idx = shift + (1:K.l);
        shift = shift + K.l;
        keep(idx) = x(idx)>0;
        K.l = sum(keep(idx));
    end
    
    % Zero out rows and columns of PSD matrices with zero diagonal entry
    if isfield(K,'s') && ~isempty(K.s)
        
        for k = 1:length(K.s)
            % Make PSD matrix symbolically
            idx = shiftPSD + ( 1 : K.s(k)^2 );
            S = reshape(x(idx),K.s(k),K.s(k));
            % Find entries to be zeroed
            d = find(diag(S)==0);
            nd = length(d);
            if nd > 0
                [i,j] = meshgrid(d,1:K.s(k));
                pos = shiftPSD + sub2ind([K.s(k),K.s(k)], [i(:); j(:)], [j(:); i(:)]);
                keep(unique(pos)) = false;
                % Update cone size
                K.s(k) = K.s(k) - nd;
            end
            % Update shift
            shiftPSD = shiftPSD + K.s(k)^2;
        end
        
        % Are we done
        if all(keep)
            done = true;
        else
            % Update At matrix and c, then project again to see if we have
            % revealed more structure
            fprintf('Removings %i variables\n',sum(~keep))
            At = At(keep,:);
            c = c(keep,:);
            [~,~,x0,PROJ,done] = removeConstraintsX(At.',b);
        end
        
    end
    
end % End of while loop

% Update some other YALMIP stuff
if isfield(K,'rank'); K.rank = K.s; end
if isfield(K,'dualrank'); K.dualrank = K.s; end
F_struc = [c, -At];

%%% END FUNCTION
end
