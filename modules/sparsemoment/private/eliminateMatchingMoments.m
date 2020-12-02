function [At,b,c,K,PROJ,bshift,y0] = eliminateMatchingMoments(At,b,c,K)

% Eliminate moments using constraints that involve at most TWO moments
% (preserve sparsity!)

% Initalize projection matrix for moments
num_moments = size(At,2);
PROJ = speye(num_moments);
bshift = 0;
y0 = sparse(num_moments,1);

% Do nothing if no equality constraints
if ~isfield(K,'f') || isempty(K.f) || K.f==0
    return
end

% Split into equality and other cones
Ateq = At(1:K.f,:); Atother = At(K.f+1:end,:);
ceq = c(1:K.f); cother = c(K.f+1:end,:);

% Find constraints that only depend on two moments
twovar = sum( spones(Ateq), 2 )<=2;
while any(twovar)
    At_twomoments = Ateq(twovar,:);
    c_twomoments = ceq(twovar);
    Ateq_others = Ateq(~twovar,:);
    ceq_others = ceq(~twovar);
    
    % Dependent constraints? Ignore them...
    indep_idx = lirows([c_twomoments, At_twomoments]);
    if ~all(indep_idx)
        c_twomoments = c_twomoments(indep_idx);
        At_twomoments = At_twomoments(indep_idx,:);
    end
    
    % Get a tolerance for cleaning
    [m,n] = size(At_twomoments);
    TOL = max(m,n)*eps(class(At_twomoments))*norm(At_twomoments,inf);

    % Now solve for the moments as y+N*z, where y is known and N satisfies 
    % At_twomoments*N = 0. We hope that N has only one entry on each row.
    y = At_twomoments\c_twomoments;
    N = fastnull((At_twomoments),'r');
    P = sum(spones(N),2);
    if ~all(P<=1)
        warning('Cannot eliminate boundary moments in a trivial way!')
        return
    end
    ceq = ceq_others - Ateq_others*y;
    Ateq = Ateq_others*N;
    Ateq = clean(Ateq, TOL);
    ceq = clean(ceq, TOL);
    
    % Eliminate zero rows?
    zeroRows = sum(spones(Ateq),2)==0;
    if any(zeroRows)
        if any(ceq(zeroRows))
            error('Inconsistent constraints: requested C = 0 with C~=0!')
        else
            Ateq = Ateq(~zeroRows,:);
            ceq = ceq(~zeroRows);
        end
    end
    
    % Update other constraints
    cother = cother - Atother*y;
    Atother = Atother*N;
    Atother = cleanzeros(Atother, TOL);
    cother = cleanzeros(cother, TOL);
    
    % Update the cost
    bshift = bshift + b.'*y;
    b = (b.'*N).';
    
    % Update projection
    y0 = y0 + PROJ*y;
    PROJ = PROJ*N;
    
    % Update number of equality constraints in this model
    K.f = length(ceq);
    
    % Assemble
    At = [Ateq; Atother];
    c = [ceq; cother];
    
    % Did we reveal anything else?
    twovar = sum( spones(Ateq), 2)<=2;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning function for sparse matrices
function A = cleanzeros(A, TOL)
    % Clean values smaller than a tolerance
    if issparse(A)
        [m,n] = size(A);
        [iA,jA,vA] = find(A);
        idx = abs(vA)>=TOL;
        A = sparse(iA(idx),jA(idx),vA(idx),m,n);
        
    else
        A(abs(A)<TOL)=0;
    end
    

end