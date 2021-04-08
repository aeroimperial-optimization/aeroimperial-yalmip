function output = callsdpagmp(interfacedata)

% CALLSDPAGMP.m Call SDPA-GMP from YALMIP
% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    23/08/2016
%
%     Copyright (C) 2016  Giovanni Fantuzzi
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ----------------------------------------------------------------------- %

% Check if max files have already been compiled
persistent iscompiled
if isempty(iscompiled); iscompiled = 0; end
if ~iscompiled
    ismex = exist('sdpagmp_read_output');
    if ismex~=3
        S = which('sdpagmp_read_output.cpp');
        compileCommand = ['mex -silent -largeArrayDims -outdir ',fileparts(S),' ',S];
        eval(compileCommand)
    end
    iscompiled = 1;
end

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% Convert from internal (sedumi) format
[mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);

if options.verbose==0
    options.sdpa_gmp.print = 'no';
else
    options.sdpa_gmp.print = 'display';
end

if options.savedebug
    ops = options.sdpa_gmp;
    save sdpa_gmpdebug mDIM nBLOCK bLOCKsTRUCT c F ops
end

if options.showprogress
    showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL SDPA-GMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solvertime = tic;
[objVal,x,X,Y,INFO] = sdpagmp(mDIM,nBLOCK,bLOCKsTRUCT,c,F,[],[],[],options.sdpa_gmp);
solvertime = toc(solvertime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here onwards, like in YALMIP native callsdpa
% Create variables in YALMIP internal format
Primal = x;

Dual = [];
for i = 1:length(Y)
    Dual = [Dual;Y{i}(:)];
end

Slack = [];
if options.saveduals
    for i = 1:length(X)
        Slack = [Slack;X{i}(:)];
    end
end

switch (INFO.phasevalue)
    case 'pdOPT'
        problem = 0;
    case {'noINFO','pFEAS','dFEAS'}
        problem = 3;
    case {'pdFEAS'}
        problem = 4;
    case 'pFEAS_dINF'
        problem = 2;
    case 'pINF_dFEAS'
        problem = 1;
    case 'pUNBD'
        problem = 2;
    case 'dUNBD'
        problem = 1;
    case 'pdINF'
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput
    solveroutput.objVal = objVal;
    solveroutput.x = x;
    solveroutput.X = X;
    solveroutput.Y = Y;
    solveroutput.INFO = INFO;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.mDIM = mDIM;
    solverinput.nBLOCK=nBLOCK;
    solverinput.bLOCKsTRUCT=bLOCKsTRUCT;
    solverinput.c=c;
    solverinput.F=F;
else
    solverinput = [];
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);

% END MAIN FUNCTION
end



% ====================================================================== %
function [objVal,x,X,Y,INFO]=sdpagmp(mDIM,nBLOCK,bLOCKsTRUCT,c,F,x0,X0,Y0,OPTION)
%
% Compute the solution of standard SDP.
% At the moment x0, X0 and Y0 are ignored as inputs. They are technically
% optional as inputs (as well as options).
%
% [objVal,x,X,Y,INFO] = sdpagmp(mDIM,nBLOCK,bLOCKsTRUCT,c,F,
%                               x0,X0,Y0,OPTION);
%
% <INPUT>
% - mDIM       : integer   ; number of primal variables
% - nBLOCK     : integer   ; number of blocks of F
% - bLOCKsTRUCT: vector    ; represetns the block structure of F
% - c          : vector    ; coefficient vector
% - F          : cell array; coefficient matrices
% - x0,X0,Y0   : cell array; initial point
% - OPTION     : structure ; options
%
% <OUTPUT>
% - objVal: [objValP objValD]; optimal value of P and D
% - x     : vector           ; optimal solution
% - X,Y   : cell arrray      ; optimal solutions
% - INFO  : structure        ; infomation of the solution
%

% SDPAGMP.m Call SDPA-GMP from YALMIP
% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    23/08/2016
%
%     Copyright (C) 2016  Giovanni Fantuzzi
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ----------------------------------------------------------------------- %
t = cputime;

if (nargin < 5 || nargin > 9)
    error('incorrect number of input arguments')
elseif nargin == 5
    % make initial points empty
    x0=[]; X0=[]; Y0=[];
    % load default parameters
    OPTION=sdpsettings;
elseif nargin == 6
    % use OPTION given by arguments
    OPTION = x0;
    % make initial points empty
    x0=[];X0=[];Y0=[];
elseif nargin == 8
    % make initial points empty
    x0=[];X0=[];Y0=[];
    % load default parameters
    OPTION=sdpsettings;
elseif nargin == 9
    % make initial points empty
    x0=[];X0=[];Y0=[];
end
ops=OPTION;

% Name of various files
inputSDPA  = [ops.inputName,'.dat-s'];
outputSDPA = [ops.outputName,'.out'];
paramSDPA = [ops.paramsName,'.sdpa'];
logSDPA = [ops.inputName,'.log'];

% If it exists, file param.sdpa will override anything written in options
% (for back compatibility), but otherwise a param.sdpa file will be created
% from the information in options and then deleted in the end.
cleanup = 0;
if ~exist([pwd,filesep,paramSDPA],'file')
    cleanup = 1;
    % write it with default parameters, otherwise failure!
    fID = fopen([pwd,filesep,paramSDPA],'w');
    fprintf(fID,'%s     unsigned int    maxIteration;           \n',int2str(ops.maxIteration));
    fprintf(fID,'%s     double          0.0 < epsilonStar;      \n',num2str(ops.epsilonStar));
    fprintf(fID,'%s     double          0.0 < lambdaStar;       \n',num2str(ops.lambdaStar));
    fprintf(fID,'%s     double          1.0 < omegaStar;        \n',num2str(ops.omegaStar));
    fprintf(fID,'%s     double          lowerBound;             \n',num2str(ops.lowerBound));
    fprintf(fID,'%s     double          upperBound;             \n',num2str(ops.upperBound));
    fprintf(fID,'%s     double          0.0 <= betaStar <  1.0; \n',num2str(ops.betaStar));
    fprintf(fID,'%s     double          0.0 <= betaBar  <  1.0, betaStar <= betaBar;\n',num2str(ops.betaBar));
    fprintf(fID,'%s     double          0.0 < gammaStar <  1.0; \n',num2str(ops.gammaStar));
    fprintf(fID,'%s     double          0.0 < epsilonDash;      \n',num2str(ops.epsilonDash));
    fprintf(fID,'%s     precision                               \n',int2str(ops.precision));
    fclose(fID);
end

% Write SDPA-GMP input file
gensdpagmpfile(inputSDPA,mDIM,nBLOCK,bLOCKsTRUCT,c,F);

% Run command in system
if ispc %For PCs via Bash on Ubuntu on Windows (called through cmd)
    pwdbash = strrep(pwd,'\','/');
    pwdbash = strrep(pwdbash,' ','\ ');
    pwdbash = strcat(replaceBetween(pwdbash,1,2,strcat('/mnt/',lower(extractBefore(pwdbash,2)))),'/');
    inputPC = [pwdbash,inputSDPA];
    outputPC = [pwdbash,outputSDPA];
    paramsPC = [pwdbash,paramSDPA];
    bashcommand = [path2sdpagmp(),'sdpa_gmp -ds ',inputPC,' -o ',outputPC,' -p ',paramsPC];
    if strcmp(ops.print,'no')
        bashcommand = [bashcommand,' >> ',pwdbash,logSDPA];
    end
    cmdcommand=['bash -c "',bashcommand,'" &'];
    %Solve SDP
    system(cmdcommand);
else %For UNIX and Macs
    shellcommand=[path2sdpagmp(),'sdpa_gmp -ds ',inputSDPA,' -o ',outputSDPA,' -p ',paramSDPA];
    if strcmp(ops.print,'no')
        %Solve SDP
        shellcommand=[shellcommand,' >> ',logSDPA];
        system(shellcommand);
    else
        %Solve SDP
        system(shellcommand);
    end
end

%To prevent weird bug before importing, it is better to pause the
%computation for 1 second. Otherwise the import could be random. This bug
%has been detected at least in Windows machines.
pause(1);

% Import result
[objVal,x,X,Y,INFO] = sdpagmp_read_output(outputSDPA,mDIM,nBLOCK,full(bLOCKsTRUCT));
pause(1);

% Clean up tmp files created in this directory
delete(inputSDPA);
delete(outputSDPA);
if cleanup
    %delete('param.sdpa');
    delete(paramSDPA);
end

INFO.cpusec = cputime-t;

% END FUNCTION
end

% ======================================================================= %
function gensdpagmpfile(filename,mDIM,nBLOCK,bLOCKsTRUCT,c,F)
%
% Generate SDP data file with SDPA Sparse format
%
% gensdpagmpfile(filename,mDIM,nBLOCK,bLOCKsTRUCT,c,F)
%
% <INPUT>
% - filename   : string    ; generated filename
% - mDIM       : integer   ; number of decision variables
% - nBLOCK     : integer   ; number of blocks of F
% - bLOCKsTRUCT: vector    ; represetns the block structure of F
% - c          : vector    ; coefficient vector
% - F          : cell array; coefficient matrices
%

% This file is a component of SDPA
% Copyright (C) 2004-2013 SDPA Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% SDPA-M: $Revision: 6.2 $
% $Id: gensdpafile.m,v 6.2 2005/05/28 02:36:40 drophead Exp $

% check the validity of arguments
if nargin ~= 6
    error('input arguments must be 6.');
elseif ~ischar(filename)
    error('1st argument must be a filename.');
end

%Giovanni Fantuzzi and Federico Fuentes on August 2017:
%YALMIP will sometimes produce these variables in sparse format so they
%should be changed to full format first for the rest of this function to
%make sense (see below calls to bLOCKsTRUCT(l), c(k) and tmpF=F{l,k};
%tmpF(i,j)).
F = cellfun(@full,F,'UniformOutput',false);
c = full(c);
bLOCKsTRUCT = full(bLOCKsTRUCT);

% open file
fid=fopen(filename, 'w');
if fid == -1
    error('Failed to open %s.', filename);
end

% comment
stamp=clock;
fprintf(fid,...
    '"Generated by gensdpagmpfile() %04d/%02d/%02d %02d:%02d:%02d"\n',...
    stamp(1),stamp(2),stamp(3),stamp(4),stamp(5),fix(stamp(6)));

% mDIM
fprintf(fid, '%d\n', mDIM);
% nBLOCK
fprintf(fid, '%d\n', nBLOCK);

% bLOCKsTRUCT
for l=1:nBLOCK
    if l~= nBLOCK
        fprintf(fid,'%d,', bLOCKsTRUCT(l));
    else
        fprintf(fid,'%d\n', bLOCKsTRUCT(l));
    end
end

% cost vector c
for k=1:mDIM
    if k ~= mDIM
        fprintf(fid,'%15.15g,', c(k));
    else
        fprintf(fid,'%15.15g\n', c(k));
    end
end

% coefficient matrices F => Sparse format
for k=1:mDIM+1
    for l=1:nBLOCK
        dim=abs(bLOCKsTRUCT(l));
        tmpF=F{l,k};
        if isempty(tmpF)
            continue;
        end
        [m,n]=size(tmpF);
        if ((m == dim && n == dim) && bLOCKsTRUCT(l) > 0)
            % normal block
            for i=1:dim
                for j=1:dim
                    if (i <= j && tmpF(i,j) ~= 0)
                        % upper triangle part only
                        fprintf(fid, '%d,%d,%d,%d,%15.15g\n',k-1,l,i,j,tmpF(i,j));
                    end
                end
            end
        elseif ((( m==1 || n==1 ) && m*n==dim) && bLOCKsTRUCT(l) < 0)
            % diagonal block with vector
            for i=1:dim
                if tmpF(i) ~= 0
                    fprintf(fid, '%d,%d,%d,%d,%15.15g\n',k-1,l,i,i,tmpF(i));
                end
            end
        elseif ((m == dim && n == dim) && bLOCKsTRUCT(l) < 0)
            % diagonal block with matrix
            for i=1:dim
                if tmpF(i,i) ~= 0
                    fprintf(fid, '%d,%d,%d,%d,%15.15g\n',k-1,l,i,i,tmpF(i,i));
                end
            end
        else
            error('Inconsistent data at F{%d,%d}',l,k);
            %fclose(fid);
            %return;
        end
    end
end

% close file
fclose(fid);

% END OF FUNCTION
end