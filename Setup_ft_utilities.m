function    varargout = Setup_ft_utilities

% SETUP_FT_UTILITIES add fieldtrip_utilities directory into path
% Usage: 
%   oldpath = SETUP_FT_UTILITIES;

% 20180517 Yuasa

nargoutchk(0,1);

rootdir = fileparts(which(mfilename));

%-- add ft_utilities
adddir = {'ft_utilities', 'fieldtrip_utilities'};
locpos = {'.','..'};
  for idir = 1:length(adddir)   % check in locpos
    addexist = false(1,length(locpos));
    for iloc = 1:length(locpos)
        addexist(idir) = exist(fullfile(rootdir,locpos{iloc},adddir{idir}),'dir');
    end
    isadd = find(addexist,1);
    if ~isempty(isadd)
        pathlist = genpath(absolute_path(fullfile(rootdir,locpos{iloc},adddir{idir})));
        break;
    end
  end
    
assert(~isempty(isadd), 'failed to find fieldtrip_utilities directory');
pathlist = strsplit(pathlist,pathsep)';

%-- remove calibration related dir
for idir = length(pathlist):-1:1
    if ~isempty(strfind(pathlist{idir},'CalibratedData')) || ~isempty(strfind(pathlist{idir},'mytransformation'))
        pathlist(idir) = [];
    end
end

%-- remove private dir
for idir = length(pathlist):-1:1
    if strfind(pathlist{idir},'ft_private')
        pathlist(idir) = [];
    end
end

%-- set low priority for ft_fix
lowPpath = {};
for idir = length(pathlist):-1:1
    if strfind(pathlist{idir},'ft_fix')
        lowPpath       = [lowPpath; pathlist(idir)];
        pathlist(idir) = [];
    end
end

%-- output
pathlist  = strjoin(pathlist,pathsep);
lowPpath  = strjoin(lowPpath,pathsep);
oldpath   = addpath(pathlist);
addpath(lowPpath,'-END');
savepath;

if nargout>0
    varargout{1} = oldpath;
end

end

function APATH = absolute_path(RPATH) 

if ~ischar(RPATH), error('Input value is not path.'); end; 

[dirname, filename, ext] = fileparts(RPATH);
[num,pathinfo] = fileattrib(dirname);
APATH = fullfile(pathinfo.Name, strcat([filename, ext]));

end