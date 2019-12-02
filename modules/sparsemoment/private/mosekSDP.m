function [x, output] = mosekSDP(A,b,c,K,options)

% A function to call mosek to solve an SDP problem in sedumi format --
% adapted from yalmip's own callmosek.m

% Get mosek format
param = options.mosek;
prob = sedumi2SDPmosek(A,b,c,K);

% Debug?
if options.savedebug
    save mosekdebug prob param
end

% Do stuff
if options.verbose == 0
    solvertime = tic;
    [r,res] = mosekopt('minimize echo(0)',prob,param);    
    solvertime = toc(solvertime);
else
    solvertime = tic;
    [r,res] = mosekopt('minimize info',prob,param);
    solvertime = toc(solvertime);
end

if res.rcode == 2001
    res.sol.itr.prosta = 'DUAL_INFEASIBLE';
end

try
    x = res.sol.itr.y;
catch
    x = nan(length(model.c),1);    
end

if options.saveduals && ~isempty(x)
    D_struc = [res.sol.itr.xx];   
    D_struc_SDP = sum(K.s.^2);
    top = 1;
    dtop = 1;
    for i = 1:length(K.s)
        X = zeros(K.s(i));
        n = K.s(i);
        I = find(tril(ones(n)));
        D_struc_SDP(dtop + I - 1) = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
        X(I) = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
        X = X + tril(X,-1)';
        D_struc = [D_struc;X(:)];
        top = top + n*(n+1)/2;
        dtop = dtop  + n^2;
    end
else
    D_struc = [];
end

switch res.sol.itr.prosta
    case 'PRIMAL_AND_DUAL_FEASIBLE'
        problem = 0;
    case 'DUAL_INFEASIBLE'
        problem = 1;
    case 'PRIMAL_INFEASIBLE'
        problem = 2;
    case 'UNKNOWN'
        try
            if isequal(res.rcodestr,'MSK_RES_TRM_STALL')
                problem = 4;
            else
                problem = 9;
            end
        catch
            problem = 9;
        end
    otherwise
        problem = -1;
end

% ---------------------------------- %
% Standard yalmip output
infostr = yalmiperror(problem,'MOSEK');
% Save solver input?
if options.savesolverinput
    solverinput.prob = prob;
    solverinput.param = options.mosek;
else
    solverinput = [];
end
% Save all data from the solver?
if options.savesolveroutput
    solveroutput.r = r;
    solveroutput.res = res;
    
else
    solveroutput = [];
end
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);

% End function
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob = sedumi2SDPmosek(A,b,c,K)

% Transform sedumi model into mosek problem -- copied from YALMIP's own
% function: yalmip2SDPmosek

prob.a = A(1:(K.f+K.l+sum(K.q)),:)';
prob.c = c(1:(K.f+K.l+sum(K.q)));

prob.blx = [-inf(K.f,1);zeros(K.l,1);-inf(sum(K.q),1)];
prob.bux = [inf(K.f+K.l+sum(K.q),1)];

prob.bardim = K.s;
prob.blc = b;
prob.buc = b;

top = 1+K.f+K.l+sum(K.q);
prob.barc.subj = [];
prob.barc.subk = [];
prob.barc.subl = [];
prob.barc.val = [];
prob.bara.subi = [];
prob.bara.subj = [];
prob.bara.subk = [];
prob.bara.subl = [];
prob.bara.val = [];

tops = [1];
for j = 1:length(K.s)
    n = K.s(j);
    tops = [tops tops(end)+n^2];
end
[ii,jj,kk] = find(A(top:top + sum(K.s.^2)-1,:));
cols = zeros(length(ii),1);
rows = zeros(length(ii),1);
allcol = [];
allrow = [];
allcon = [];
allvar = [];
allval = [];
for j = 1:length(K.s)    
    ind = find(ii>=tops(j) & ii<=tops(j+1)-1);
    iilocal = ii(ind)-tops(j)+1;
    col = ceil(iilocal/K.s(j));
    row = iilocal - (col-1)*K.s(j);
    allcol = [allcol col(:)'];
    allrow = [allrow row(:)'];
    allvar = [allvar jj(ind(:))'];
    allval = [allval kk(ind(:))'];
    allcon = [allcon repmat(j,1,length(col))];
end
keep = find(allrow >= allcol);
allcol = allcol(keep);
allrow = allrow(keep);
allcon = allcon(keep);
allvar = allvar(keep);
allval = allval(keep);
%allvar = jj(keep);
%allval = kk(keep);
prob.bara.subi = [prob.bara.subi allvar];
prob.bara.subj = [prob.bara.subj allcon];
prob.bara.subk = [prob.bara.subk allrow];
prob.bara.subl = [prob.bara.subl allcol];
prob.bara.val = [prob.bara.val allval];

for j = 1:length(K.s)
    n = K.s(j);
    Ci = c(top:top+n^2-1);
    Ci = tril(reshape(Ci,n,n));
    [k,l,val] = find(Ci);
    prob.barc.subj = [prob.barc.subj j*ones(1,length(k))];
    prob.barc.subk = [prob.barc.subk k(:)'];
    prob.barc.subl = [prob.barc.subl l(:)'];
    prob.barc.val = [prob.barc.val val(:)'];     
    top = top + n^2;
end

if K.q(1)>0    
    prob.cones.type   = [repmat(0,1,length(K.q))];
    top = 1 + K.f + K.l;
    prob.cones.sub = [];
    prob.cones.subptr = [];
    for i = 1:length(K.q)
        prob.cones.subptr = [ prob.cones.subptr 1+length(prob.cones.sub)];
        prob.cones.sub = [prob.cones.sub top:top+K.q(i)-1];
        top = top + K.q(i);
    end
end

% End function
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%