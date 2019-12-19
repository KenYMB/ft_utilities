function [cfg] = ft_connectivityplot_twocondition(cfg, varargin)

% FT_CONNECTIVITYPLOT plots channel-level frequency resolved connectivity. The
% data are rendered in a square grid of subplots, each subplot containing the
% connectivity spectrum between the two respective channels.
%
% Use as
%   ft_connectivityplot(cfg, data)
%
% The input data is a structure containing the output to FT_CONNECTIVITYANALYSIS
% using a frequency domain metric of connectivity. Consequently the input
% data should have a dimord of 'chan_chan_freq'.
%
% The cfg can have the following options:
%   cfg.parameters  = 2xN cell-strings, the functional parameter to be plotted in upper-right & lower-left field (default = {'cohspctrm'; 'cohspctrm'})
%   cfg.titles      = 2x1 cell-strings, the figure title of each upper-right & lower-left field (default = 'auto')
%   cfg.legend      = 'yes' or 'no', if show legend (default = 'yes')
%   cfg.xytitle     = 'yes' or 'no', if show x,y titles (default = 'yes')
%   cfg.axis        = 'yes' or 'no', if show axis labels (default = 'yes')
%   cfg.xlim        = selection boundaries over first dimension in data (e.g., freq)
%                     'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim        = plotting limits for color dimension, 'maxmin', 'maxabs' or [zmin zmax] (default = 'maxmin')
%   cfg.channel     = list of channels to be included for the plotting (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.graphcolor  = string, the color list of plotted lines
%
% See also FT_CONNECTIVITYANALYSIS, FT_CONNECTIVITYSIMULATION, FT_MULTIPLOTCC, FT_TOPOPLOTCC

% Copyright (C) 2011-2013, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% using: fieldtrip, ft_private(rollback_provenance)

% 20170224 Yuasa

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

ft_private('silent');
try

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'zparam', 'parameters'});
cfg = ft_checkconfig(cfg, 'renamed', {'color',  'graphcolor'}); % to make it consistent with ft_singleplotER

% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel',   'all');
cfg.parameters = ft_getopt(cfg, 'parameters', {'cohspctrm';'cohspctrm'});
cfg.titles     = ft_getopt(cfg, 'titles', 'auto');
cfg.legend     = ft_getopt(cfg, 'legend', 'yes');
cfg.xytitle    = ft_getopt(cfg, 'xytitle', 'yes');
cfg.axis       = ft_getopt(cfg, 'axis', 'yes');
cfg.zlim       = ft_getopt(cfg, 'zlim',      'maxmin');
cfg.xlim       = ft_getopt(cfg, 'xlim',      'maxmin');
cfg.graphcolor = ft_getopt(cfg, 'graphcolor', 'brgkywbrgkywbrgkywbrgkywbrgkywbrgkywbrgkyw');

% check if the input data is valid for this function
% ensure that the input is correct
Ndata = numel(varargin);
dtype = cell(Ndata, 1);
iname = cell(Ndata+1, 1);
for k = 1:Ndata
  % check if the input data is valid for this function
  varargin{k} = ft_checkdata(varargin{k}, 'datatype', {'timelock', 'freq'});
  dtype{k}    = ft_datatype(varargin{k});

  % convert into the the supported dimord
  if strcmp(varargin{k}.dimord, 'chan_chan_freq')
    % that's ok
  elseif strcmp(varargin{k}.dimord, 'chancmb_freq')
    % convert into 'chan_chan_freq'
    varargin{k} = ft_checkdata(varargin{k}, 'cmbrepresentation', 'full');
  else
    error('the data should have a dimord of %s or %s', 'chan_chan_freq', 'chancmb_freq');
  end
  
  % this is needed for correct treatment of graphcolor later on
  if nargin>1,
    if ~isempty(inputname(k+1))
      iname{k+1} = inputname(k+1);
    else
      iname{k+1} = ['input',num2str(k,'%02d')];
    end
  else
    % not yet supported
    iname{k+1} = cfg.inputfile{k};
  end
end

if Ndata >1,
  if ~all(strcmp(dtype{1}, dtype))
    error('input data are of different type; this is not supported');
  end
end

% check cfg.parameters
if ~iscellstr(cfg.parameters) || numel(cfg.parameters) ~= 2
  if size(cfg.parameters,1) ~= 2 || size(cfg.parameters,2) ~= Ndata
    error('cfg.parametes must be 2xN cell array of strings');
  end
end
if strcmp(cfg.titles,'auto')
  if numel(cfg.parameters) == 2
    cfg.titles        = cfg.parameters;
  else
    cfg.titles        = {'';''};
  end
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'channel'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

% Get physical min/max range of x:
if ischar(cfg.xlim) && strcmp(cfg.xlim,'maxmin')
  xmin = inf;
  xmax = -inf;
  for k = 1:Ndata
    xmin = min(xmin,varargin{k}.freq(1));
    xmax = max(xmax,varargin{k}.freq(end));
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end
cfg.xlim = [xmin xmax];

% Get physical min/max range of z:
if ischar(cfg.zlim) && strcmp(cfg.zlim,'maxmin')
  zmin = inf;
  zmax = -inf;
  for k = 1:Ndata
    zmin = min(zmin,min(varargin{k}.(cfg.parameters{1})(:)));
    zmax = max(zmax,max(varargin{k}.(cfg.parameters{1})(:)));
    zmin = min(zmin,min(varargin{k}.(cfg.parameters{2})(:)));
    zmax = max(zmax,max(varargin{k}.(cfg.parameters{2})(:)));
  end
elseif ischar(cfg.zlim) && strcmp(cfg.zlim,'maxabs')
  zmax = -inf;
  for k = 1:Ndata
    zmax = max(zmax,max(abs(varargin{k}.(cfg.parameters{1})(:))));
    zmax = max(zmax,max(abs(varargin{k}.(cfg.parameters{2})(:))));
  end
  zmin = -zmax;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end
cfg.zlim = [zmin zmax];

% make the function recursive if Ndata>1
if Ndata>1
  data = varargin{1};
  tmpcfg = cfg;
  if numel(cfg.parameters)==2
    % do nothing
  else
    tmpcfg.parameters = cfg.parameters{:,1};
  end
  ft_connectivityplot_twocondition(tmpcfg, data);
  tmpcfg = cfg;
  tmpcfg.legend    = 'no';
   
  if ischar(cfg.graphcolor),        colorLabels = [iname{2} '_feedforward = ' tmpcfg.graphcolor(1) '\n' ...
                                                   iname{2} '_feedback    = ' tmpcfg.graphcolor(2) '\n'];
  elseif isnumeric(cfg.graphcolor), colorLabels = [iname{2} '_feedforward = ' num2str(tmpcfg.graphcolor(1, :)) '\n' ...
                                                   iname{2} '_feedback    = ' num2str(tmpcfg.graphcolor(2, :)) '\n'];
  end
    
  for k = 2:Ndata
    tmpcfg.xytitle   = 'no';
    tmpcfg.axis      = 'no';
  
    if ischar(cfg.graphcolor),     tmpcfg.graphcolor = tmpcfg.graphcolor(2:end);
    else isnumeric(cfg.graphcolor),tmpcfg.graphcolor = tmpcfg.graphcolor(2:end,:);
    end
    
    tmpcfg.holdfig = 1;
    if numel(cfg.parameters)==2
      % do nothing
    else
      tmpcfg.parameters = cfg.parameters{:,k};
    end
    ft_connectivityplot_twocondition(tmpcfg, varargin{k});
    
    if ischar(cfg.graphcolor);        colorLabels = [colorLabels iname{k+1} '=' tmpcfg.graphcolor(1) '\n'];
    elseif isnumeric(cfg.graphcolor); colorLabels = [colorLabels iname{k+1} '=' num2str(tmpcfg.graphcolor(1, :)) '\n'];
    end
  end
  
  if strcmp(cfg.legend, 'yes')
    ft_plot_text(0.5, (numel(varargin{k}.label)+1).*1.2-0.5, sprintf(colorLabels), 'interpreter', 'none', 'horizontalalignment', 'right');
  end

  return;
else
  data = varargin{1};
end

if ~isfield(data, cfg.parameters{1})
  error('the data does not contain the requested parameter %s', cfg.parameters{1});
end
if ~isfield(data, cfg.parameters{2})
  error('the data does not contain the requested parameter %s', cfg.parameters{2});
end

% get the selection of the data
tmpcfg           = [];
tmpcfg.frequency = cfg.xlim;
data             = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

nchan = numel(data.label);
nfreq = numel(data.freq);

% merge two conditions
data1 = rmfield(data,cfg.parameters);
data2 = data1;
data1.newspctrm     = data.(cfg.parameters{1});
data2.newspctrm     = permute(data1.newspctrm, [2 1 3]);    % flip feedforward & feedback
for k = 1:nchan
  for m = 1:nchan
    if k > m
        data1.newspctrm(k,m,:)  = data.(cfg.parameters{2})(m,k,:);      % cond1-fb <= cond2-ff
        data2.newspctrm(k,m,:)  = data.(cfg.parameters{2})(k,m,:);      % cond1-ff <= cond2-fb
    else
        % data1: remain cond1-ff
        % data2: remain cond1-fb
    end
  end
end

% call ft_connectivityplot
tmpcfg           = cfg;
tmpcfg.parameter = 'newspctrm';
tmpcfg.legend    = 'no';
tmpcfg.xytitle   = 'no';
ft_connectivityplot(tmpcfg, data1, data2);

if strcmp(cfg.legend, 'yes')
  if ischar(cfg.graphcolor),        colorLabels = ['feedforward = ' tmpcfg.graphcolor(1) '\n' ...
                                                   'feedback    = ' tmpcfg.graphcolor(2) '\n'];
  elseif isnumeric(cfg.graphcolor), colorLabels = ['feedforward = ' num2str(tmpcfg.graphcolor(1, :)) '\n' ...
                                                   'feedback    = ' num2str(tmpcfg.graphcolor(2, :)) '\n'];
  end
  ft_plot_text(0.5, (numel(varargin{1}.label)+1).*1.2-0.5, sprintf(colorLabels), 'interpreter', 'none', 'horizontalalignment', 'right');
end

% add 'from' and 'to' labels
if strcmp(cfg.xytitle, 'yes')
ft_plot_text(-0.5,           (nchan + 1)/1.7, cfg.titles{2}, 'rotation', 90);
ft_plot_text((nchan + 1)/1.7, (nchan + 1)*1.2+0.4, cfg.titles{1});
end

set(gcf, 'color', [1 1 1]);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance

ft_private('remove','silent');

catch ME
    ft_private('remove','silent');
    rethrow(ME)
end % try
