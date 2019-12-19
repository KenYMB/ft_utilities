function varargout = ft_private(option,issilent)

% oldpath = FT_PRIVATE(option)
% ft_privateフォルダへのパスを設定する
% 
% option = 
%   'set(default)'  : パスを追加
%   'remove'        : パスを削除
%   'reset'         : 直前のft_private('set')実行時の状態に戻す
% 
% FT_PRIVATE(option, 'silent')
% メッセージを表示させない

% 20160812 Yuasa
% 20160906 Yuasa: expand the load locations
% 20161221 Yuasa: enable to refer 'private' (it does not do anything)
% 20161221 Yuasa: minor fix
% 20170203 Yuasa: add silent mode & disable to refer 'private'
%                 because there might be a private directory of another aim
% 20170510 Yuasa: minor change
% 20170713 Yuasa: return current path
% 20170826 Yuasa: add 'reset' option

persistent prevprivate privatepath
if isempty(prevprivate)
    prevprivate = false;
end

if (nargin < 1) || isempty(option)
    option   = 'set';
end
if strcmpi(option,'silent')
    %-- enable usage of ft_private('silent',option)
    if nargin > 1
        option   = issilent;
    else
        option   = 'set';
    end
    issilent = true;
elseif (nargin > 1) && strcmpi(issilent,'silent')
    issilent = true;
else
    issilent = false;
end

%-- search private path
if ~isempty(privatepath) || ~exist(privatepath,'dir')
    ST = dbstack('-completenames');
    %-- check same directory
    privatepath = fileparts(getfield(ST,{1},'file'));
    privatepath = fullfile(privatepath,'ft_private');
    
    %-- check user directory
    upathlist = strsplit(userpath,pathsep);
    for ilp = 1:length(upathlist)
        if exist(privatepath,'dir'),    break;      end
        privatepath = fullfile(upathlist{ilp}, 'toolbox','ft_private');
        if exist(privatepath,'dir'),    break;      end
        privatepath = fullfile(upathlist{ilp}, 'ft_private');
    end
    if ~exist(privatepath,'dir')
        privatepath = '';
        %-- check if the function call has private directory
        if length(ST) > 1
            privateforfunc  = fullfile(fileparts(getfield(ST,{2},'file')),'private');
            if exist(privateforfunc,'dir'),     return;         end
        end
    end
    
    %-- skip error if private directory exist
    assert(~isempty(privatepath), sprintf('Please prepare ''ft_private'' folder at ''%s''.',fileparts(privatepath)));
end

%-- set path
oldpath = path;
if exist(privatepath,'dir')
    curwarn = warning;
    if issilent,    warning('off','all');   end
    switch option
        case 'set',         ws = warning('backtrace', 'off');
                            warning('Add ''ft_private'' to search path.');
                            warning(ws); % return to the previous warning level
                            addpath(privatepath);
                            if match_str(strsplit(oldpath,pathsep),privatepath)
                             prevprivate = true;
                            else
                             prevprivate = false;
                            end
                            
        case 'remove',      ws = warning('backtrace', 'off');
                            warning('Remove ''ft_private'' from search path.');
                            warning(ws); % return to the previous warning level
                            rmpath(privatepath);
                            prevprivate = false;
        case 'reset',
          if ~prevprivate,  ws = warning('backtrace', 'off');
                            warning('Remove ''ft_private'' from search path.');
                            warning(ws); % return to the previous warning level
                            rmpath(privatepath);
                            prevprivate = false;
          end

    end
    warning(curwarn);
end
if nargout>0,   varargout{1} = oldpath; end
