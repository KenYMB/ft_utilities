function ft_fix(varargin)

% FT_FIX
%   フォルダへのパスを最優先に設定
%   ft_fix('end') で優先度最低に設定
% 
% See also, FT_DEFAULTS

% 2014 Yuasa: wrote
% 20170216 Yuasa: enable './ft_fix'
% 20170510 Yuasa: minor change
% 20170626 Yuasa: add 'end' option
% 20170725 Yuasa: add ft_fix_flag for the case it is called by functions

persistent ft_fix_flag
if isempty(ft_fix_flag), ft_fix_flag=0; end     % set default
if length(dbstack)<=1,   ft_fix_flag=0; end     % reset when it is called manually

func_path = fileparts(mfilename('fullpath'));
    % save as "getfield(dbstack('-completenames'),{1},'file')"

if exist(fullfile(func_path,'ft_fix'),'dir')
    path_dir  = fullfile(func_path,'ft_fix');
else
    %-- set the including directory as the ft_fix path
    path_dir    = func_path;
end

if nargin>0 && (strcmpi(varargin{1},'end') || strcmpi(varargin{1},'-end'))
    %-- set ft_fix at the bottom of the search path
    if ft_fix_flag<=1
      addpath(path_dir,'-end');
    end
    ft_fix_flag=max(ft_fix_flag-1,0);
else
    %-- set ft_fix at the top of the search path
    addpath(path_dir);
    ft_fix_flag=ft_fix_flag+1;
end