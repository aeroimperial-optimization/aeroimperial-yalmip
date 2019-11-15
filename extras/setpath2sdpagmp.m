function setpath2sdpagmp(varargin)

% setpath2sdpagmp
%
% Set the path to the sdpa_gmp executable.
% With no input arguments, it is assumed that the executable sdpa_gmp
% is installed in the location
%
% /usr/local/bin/
%
% Otherwise, you can set the path to the executable by calling
%
% >> setpath2sdpagmp(PATH2EXE)
%
% where PATH2EXE is a string specifying the path to the sdpa_gmp executable.


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


% ----------------------------------------------------------------------- %
% CHECK INPUTS
% ----------------------------------------------------------------------- %
if nargin == 0
    PATH2EXE = fileparts('/usr/local/bin/');%Assumes sdpa_gmp in /usr/local/bin/
elseif nargin == 1 && ischar(varargin{1})
    PATH2EXE = fileparts(varargin{1});
elseif nargin == 1 && ~ischar(varargin{1})
    error('Input must be a string specifying the path to the sdpa_gmp executable.')
else
    error('Too many inputs!')
end

% Check it is actually the correct location
if ispc
    pwdbash=strrep(pwd,'\','/');
    pwdbash=strrep(pwdbash,' ','\ ');
    pwdbash=strcat(replaceBetween(pwdbash,1,2,strcat('/mnt/',lower(extractBefore(pwdbash,2)))),'/');
    bashcommand=['[ -f ',PATH2EXE,'/sdpa_gmp ] && touch ',pwdbash,'1.txt || :'];
    cmdcommand=['bash -c "' bashcommand '" &'];
    system(cmdcommand);
    existfile=exist('1.txt','file');
    if existfile==2
        delete('1.txt');
    end
else
    existfile=exist([PATH2EXE,filesep,'sdpa_gmp'],'file');
end
if existfile~=2
    str=['Cannot find the sdpa_gmp executable at the directory ''',PATH2EXE,'/''. '];
    str=[str,'Please make sure you call install_sdpa_gmp(''/path/to/sdpa/gmp/''), where '];
    str=[str,'/path/to/sdpa/gmp/ is a UNIX-based path to either the sdpa_gmp executable '];
    str=[str,'itself or to the directory, ending with ''/'', where the sdpa_gmp executable lies.'];
    error(str);
end

% Update path in sdpsettings.m
fname = fileparts(which('sdpsettings'));
fname = [fname,filesep,'path2sdpagmp.m'];
A = regexp( fileread(fname), '\n', 'split');

% Add the information to sdpsettings.m
for i = 1:length(A)
    substr='path ='; 
    if ~isempty(regexp(A{i},substr,'once'))
        A{i} = sprintf('path = ''%s'' ;',[PATH2EXE,'/']);
        break
    end
end

fid = fopen(fname, 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);

% ----------------------------------------------------------------------- %
% CLEAR YALMIP CACHED SOLVERS
% ----------------------------------------------------------------------- %
clear('compileinterfacedata.m')

end