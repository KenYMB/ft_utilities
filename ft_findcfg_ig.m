function [val, status] = ft_findcfg_ig(cfg, var, ignore, ignoreempty)
% FT_FINDCFG_IG searches for an element in the cfg structure
% or in the nested previous cfgs
%
% Use as
%   val = ft_findcfg(cfg, var, {ignore1, ignore2})
% where the name of the variable should be specified as string.
% if 'val' is equal to one of 'ingnores', ft_findcfg continues to search
% 
%   val = ft_findcfg(cfg, var, [], 1)
% if fourth input argument is true, ft_findcfg continues to search for
% empty value
%
% e.g.
%   trl   = ft_findcfg(cfg, 'trl')
%   event = ft_findcfg(cfg, 'event')
%
% See also FT_FINDCFG

% Copyright (C) 2006, Robert Oostenveld

% 20170817 Yuasa

val   = [];
depth = 0;
status = 0;

if nargin < 3
    ignore = {};
end
if ~iscell(ignore)
    ignore = {ignore};
end

if nargin < 4
    ignoreempty = false;
else
    ignoreempty = logical(ignoreempty(1));
end

while ~status
  depth = depth + 1;
  if ~isempty(cfg)
    if issubfield(cfg,  var)
      val = getsubfield(cfg, var);
      status = 1;
      for iig = 1:numel(ignore)
        try
          if (~isempty(val) && isequal(val,ignore{iig})) ...
                || (ignoreempty && isempty(val))
            val    = [];
            status = 0;
          end
        end
      end
    elseif issubfield(cfg, 'previous');
      [val, status] = ft_findcfg_ig(cfg.previous, var, ignore, ignoreempty);
      if status, break; end;
    elseif iscell(cfg)
      for i=1:length(cfg)
        [val, status] = ft_findcfg_ig(cfg{i}, var, ignore, ignoreempty);
        if status, break; end;
      end
    else
      status = -1;
      break
    end
  else
    status = -1;
    break
  end
end
