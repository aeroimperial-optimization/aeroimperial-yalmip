function [C,D] = corrsparsity(exponent_p_monoms,options)

% Detect correlative sparsity of a polynomial with specified support. Use
% the chordal extension approach if needed.
if options.sos.csp
    
    % Display?
    if options.verbose > 0
        fprintf('Detecting correlative sparsity...')
    end
    
    % Build CSP matrix
    n = size(exponent_p_monoms,2);
    C = zeros(n,n);
    for i = 1:size(exponent_p_monoms,1)
        row = find(exponent_p_monoms(i,:));
        for j = 1:length(row)
            C(row(j),row(j)) = 1;
            for k = 2:length(row)
                C(row(j),row(k)) = 1;
                C(row(k),row(j)) = 1;
            end
        end
    end
    
    % Extend? If so, ensure C is diagonally dominant so Cholesky exists.
    % When Cholesky-ing, use the specified reordering
    if options.sos.csp_chordal_extension
        Cnorm = norm(C,1);
        C = chol(C+Cnorm*eye(n));
        if strcmpi(options.sos.csp_ordering,'amd')
            s = symamd(C);
        elseif strcmpi(options.sos.csp_ordering,'rcm')
            s = symrcm(C);
        else
            s = 1:n;
        end
        C = chol(C(s,s));
    end 
    
    % Detect sparsity
    for i = 1:size(C,1)
        col = s(C(i,:)>0);
        if i>1
            is_in = 0;
            for j = 1:length(D)
                if all(ismember(col,D{j}))
                    is_in = 1;
                end
            end
            if ~is_in
                D{end+1} = col;
            end
        else
            D{1} = col;
        end
    end
    if length(D)>1 && options.verbose>0
        [uu,~,oo] = unique(cellfun('prodofsize',D));
        for i = 1:length(uu)
            n_this = length(find(oo==i));
            the_text = [num2str(uu(i)) '(' num2str(n_this) ')' ' '];
        end
        fprintf([the_text,'\n']);
    end
else
    C = ones(size(exponent_p_monoms,2));
    D{1} = 1:size(exponent_p_monoms,2);
end
