function    varargout = Setup_ft_utilities

% SETUP_FT_UTILITIES add fieldtrip_utilities directory into path
% Usage: 
%   oldpath = SETUP_FT_UTILITIES;

% 20180517 Yuasa
% 20251031 Yuasa: major update on the recent environment

nargoutchk(0,1);

%-- add ft_utilities
rootdir = fileparts(which(mfilename));
pathlist = genpath(rootdir);
pathlist = strsplit(pathlist,pathsep)';
filesepexp =  regexptranslate('escape',filesep);

%-- remove hidden dir
idx = ~cellfun(@isempty,regexp(pathlist,[filesepexp '\..'],'once'));
pathlist(idx) = [];

%-- remove private dir
idir = ~cellfun(@isempty,regexp(pathlist,[filesepexp 'ft_private'],'once'));
pathlist(idir) = [];

%-- set low priority for ft_fix
idir = ~cellfun(@isempty,regexp(pathlist,[filesepexp 'ft_fix'],'once'));
lowPpath = pathlist(idir);
pathlist(idir) = [];

%-- output
pathlist  = strjoin(pathlist,pathsep);
lowPpath  = strjoin(lowPpath,pathsep);
oldpath   = addpath(pathlist);
addpath(lowPpath,'-END');
% savepath;

if nargout>0
    varargout{1} = oldpath;
end

end