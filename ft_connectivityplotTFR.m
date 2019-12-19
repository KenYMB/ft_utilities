function [cfg] = ft_connectivityplotTFR(cfg, varargin)

% FT_CONNECTIVITYPLOT3D plots channel-level frequency resolved connectivity. The
% data are rendered in a square grid of subplots, each subplot containing the
% connectivity spectrum between the two respective channels.
%
% Use as
%   ft_connectivityplotTFR(cfg, data)
%
% The input data is a structure containing the output to FT_CONNECTIVITYANALYSIS
% using a frequency domain metric of connectivity. Consequently the input
% data should have a dimord of 'chan_chan_freq_time'.
%
% The cfg can have the following options:
%   cfg.parameter   = string, the functional parameter to be plotted (default = 'cohspctrm')
%   cfg.xlim        = selection boundaries over first dimension in data (e.g., time)
%                     'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim        = selection boundaries over first dimension in data (e.g., freq)
%                     'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim        = plotting limits for color dimension, 'maxmin', 'maxabs' or [zmin zmax] (default = 'maxmin')
%   cfg.colorbar    = 'yes', 'no' (default = 'yes')
%   cfg.colormap    = any sized colormap or colormap name (default = 'jet')
%   cfg.channel     = list of channels to be included for the plotting (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.fontsize    = font size of this figure (default = 10)
% 
% The following options are also available
%   cfg.maskparameter  = field in the data to be used for masking of data
%                        (not possible for mean over multiple channels, or when input contains multiple subjects
%                        or trials)
%   cfg.maskstyle      = style used to masking, 'opacity', 'saturation' or 'outline' (default = 'opacity')
%                        use 'saturation' or 'outline' when saving to vector-format (like *.eps) to avoid all sorts of image-problems
%   cfg.maskalpha      = alpha value between 0 (transparant) and 1 (opaque) used for masking areas dictated by cfg.maskparameter (default = 1)
% 
% The following option is the switch to call FT_CONNECTIVITYPLOT to draw
% time-domain plot
%   cfg.toi         = vector, call FT_CONNECTIVITYPLOT for each toi
% 
% The following is the additional option for FT_CONNECTIVITYPLOT to plot
% multiple data simultaneously
%   cfg.legend      = 1xN cell-array, legend for each plot
%
% See also FT_CONNECTIVITYPLOT, FT_CONNECTIVITYANALYSIS, FT_CONNECTIVITYSIMULATION, FT_MULTIPLOTCC, FT_TOPOPLOTCC

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

% 20161220 Yuasa: create based on ft_connectivityplot
% 20170117 Yuasa: minor bug fix
% 20170209 Yuasa: adjust axis position
% 20170328 Yuasa: suppress nargchk warning in uimage, enable cfg.toi
% 20170403 Yuasa: minor fix
% 20170430 Yuasa: enable mask
% 20171204 Yuasa: minor fix
% 20180406 Yuasa: enable cfg.legend
% 20180725 Yuasa: minor update

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

% set private functions
ft_private('set','silent')
% set off nargchk warning for uimage
curwarn = warning('off','MATLAB:nargchk:deprecated');
try

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamed', {'color',  'graphcolor'}); % to make it consistent with ft_singleplotER

% set the defaults
cfg.channel   = ft_getopt(cfg, 'channel',   'all');
cfg.parameter = ft_getopt(cfg, 'parameter', 'cohspctrm');
cfg.zlim      = ft_getopt(cfg, 'zlim',      'maxmin');
cfg.xlim      = ft_getopt(cfg, 'xlim',      'maxmin');
cfg.ylim      = ft_getopt(cfg, 'ylim',      'maxmin');
% cfg.graphcolor = ft_getopt(cfg, 'graphcolor', 'brgkywrgbkywrgbkywrgbkyw');    % change to follow colormap
cfg.colorbar   = ft_getopt(cfg, 'colorbar',     'yes');
cfg.colormap   = ft_getopt(cfg, 'colormap',     'jet');
cfg.maskalpha      = ft_getopt(cfg, 'maskalpha',     1);
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter', []);
cfg.maskstyle      = ft_getopt(cfg, 'maskstyle',    'opacity');
cfg.masknans       = ft_getopt(cfg, 'masknans',     'yes');
cfg.fontsize       = ft_getopt(cfg, 'fontsize',     10);
cfg.toi            = ft_getopt(cfg, 'toi', []);

Ndata       = numel(varargin);
haslegend   = isfield(cfg,'legend') && length(cfg.legend) >= Ndata;

if ischar(cfg.colormap)
    htmp = figure('Visible','off');
    cfg.colormap = colormap(cfg.colormap);
    delete(htmp);
end
cfg.graphcolor = ft_getopt(cfg, 'graphcolor', cfg.colormap(round(linspace(1,length(cfg.colormap),Ndata)),:));


% check if the input data is valid for this function
% ensure that the input is correct
dtype = cell(Ndata, 1);
iname = cell(Ndata+1, 1);
for k = 1:1  % skip the process for Ndata>1
  % check if the input data is valid for this function
  varargin{k} = ft_checkdata(varargin{k}, 'datatype', {'timelock', 'freq'});
  dtype{k}    = ft_datatype(varargin{k});
  dimord      = getdimord(varargin{k}, cfg.parameter);
  dimtok      = tokenize(dimord, '_');
  if length(dimtok) ~= ndims(varargin{k}.(cfg.parameter))
      dimord  = varargin{k}.dimord;
  end

  % convert into the the supported dimord
  dim_chan           = find(strcmp(dimtok,'chan'));
  dim_chancmb        = find(strcmp(dimtok,'chancmb'));
  dim_freq           = find(strcmp(dimtok,'freq'));
  dim_time           = find(strcmp(dimtok,'time'));
    
  if length(dim_chan)==2 && length(dim_freq)==1 && length(dim_time)==1
    % that's ok
    varargin{k}.(cfg.parameter)   = permute(varargin{k}.(cfg.parameter), [dim_chan dim_chancmb min(dim_freq,dim_time) max(dim_freq,dim_time)]);
  elseif length(dim_chancmb)==1 && length(dim_freq)==1 && length(dim_time)==1
    % convert into 'chan_chan_freq_time'
    varargin{k}.(cfg.parameter)   = permute(varargin{k}.(cfg.parameter), [dim_chan dim_chancmb min(dim_freq,dim_time) max(dim_freq,dim_time)]);
    varargin{k} = ft_checkdata(varargin{k}, 'cmbrepresentation', 'full');
  elseif strcmp(varargin{k}.dimord, 'chan_chan_freq') || strcmp(varargin{k}.dimord, 'chancmb_freq')
    warning(sprintf('the data doesn''t have a dimension of time\n\t  forwarding the data to ft_connectivityplot ...'))
    ft_private('remove','silent')
    ft_fix;
    try
    tmpcfg = cfg;
    if ~strcmpi(tmpcfg.ylim,'maxmin')      % forwarding cfg.ylim to cfg.xlim (consider cfg.ylim is specified for freq)
        tmpcfg.xlim = tmpcfg.ylim;
    end
    tmpcfg.legend = 'yes';
    ft_connectivityplot(tmpcfg, varargin{:});
    catch ex
    ft_fix('end');
    rethrow(ex);
    end
    ft_fix('end');
    legend('off'); colorbar('off');
    % set plot name
    if nargin > 2
      hf = gcf;
      texlist  = hf.Children.Children.findobj('Type','text');
      texstrng = {texlist.String};
      islegend = zeros(length(texstrng),1);
      for ilp = 1:length(texstrng)
        islegend(ilp) = iscellstr(texstrng{ilp}) && min(strfind([texstrng{ilp}{:}],'input'));
      end
      hleg = texlist(find(islegend,1));
        linH = hleg.Extent(4)/(Ndata+1);
        linX = hleg.Position(1)+[-0.25 0];
        if ~ischar(cfg.graphcolor)
          set(get(hleg,'Parent'),'Clipping','off');
          set(hleg,'Clipping','off','Interpreter','tex');
        end
      for ilp = 1:Ndata
          linY = hleg.Extent(2) + (Ndata-ilp+1.5)*linH;
        if haslegend
           if isempty(cfg.legend{ilp})
             hleg.String{ilp} = '';
           elseif ischar(cfg.graphcolor)
             hleg.String{ilp} = sprintf('%s%s',cfg.legend{ilp},hleg.String{ilp}(8:end));
           else
             hleg.String{ilp} = sprintf('\\color{black} %s   \\color[rgb]{%f %f %f}%s%s%s%s',...
                                        cfg.legend{ilp},cfg.graphcolor(ilp,:),[8212 8212 8212 8212]);
%              hleg.String{ilp} = sprintf('%s           ',cfg.legend{ilp});
%              hl=ft_plot_line( linX, [linY linY], 'color', cfg.graphcolor(ilp,:),'Clipping','off');
           end
        elseif ~isempty(inputname(ilp+1))
           if ischar(cfg.graphcolor)
             hleg.String{ilp} = sprintf('%s%s',inputname(ilp+1),hleg.String{ilp}(8:end));
           else
             hleg.String{ilp} = sprintf('\\color{black} %s   \\color[rgb]{%f %f %f}%s%s%s%s',...
                                        cfg.legend{ilp},cfg.graphcolor(ilp,:),[8212 8212 8212 8212]);
%              hleg.String{ilp} = sprintf('%s           ',inputname(ilp+1));
%              hl=ft_plot_line( linX, [linY linY], 'color', cfg.graphcolor(ilp,:),'Clipping','off');
           end
        end
      end
    end
    return
  else
    error('the data should have a dimord of %s or %s', 'chan_chan_freq_time', 'chancmb_freq_time');
  end
  isfull = true;
  haslabelcmb = false;
  
  % set x_y param
  xparam = dimtok{max(dim_freq,dim_time)};
  yparam = dimtok{min(dim_freq,dim_time)};
  
  % ft_connectivityplot for toi
  if ~isempty(cfg.toi)
    % skip the process for Ndata>1
    if nargin > 2, warning('only first data is plotted.');  end
    tempdata = rmfield(varargin{1},'time');
    tempdata.dimord = 'chan_chan_freq';
    tempdata.(cfg.parameter) = permute(tempdata.(cfg.parameter),[1 2 dim_freq dim_time]);
    freqdata  = cell(1,length(cfg.toi));
    for tlp = 1:length(cfg.toi)
        tidx = nearest(varargin{1}.time,cfg.toi(tlp));
        freqdata{tlp} = tempdata;
        freqdata{tlp}.(cfg.parameter) = tempdata.(cfg.parameter)(:,:,:,tidx);
    end
    ft_private('remove','silent')
    ft_fix;
    try
    tmpcfg = cfg;
    tmpcfg.legend = 'yes';
    ft_connectivityplot(tmpcfg, freqdata{:});
    catch ex
    ft_fix('end');
    rethrow(ex);
    end
    ft_fix('end');
    legend('off'); colorbar('off');
    % set plot name
    wint  = max(1,fix(log10(nanmax(cfg.toi(:)))+1));
    wdeci = max(0,size(num2str(reshape(cfg.toi-fix(cfg.toi),[],1)),2)-2);  % -2: subtract '0.'
      fleg = min(6,wint+logical(wdeci)+wdeci+1);          % logical(wdeci): '.', +1:' ', [1 6]
      wleg = max(0,min(fleg-wint-2,wdeci));               % -2: ' '&'.', [0 4]
     if numel(cfg.toi)>1
         hf = gcf;
         texlist  = hf.Children.Children.findobj('Type','text');
         texstrng = {texlist.String};
         islegend = zeros(length(texstrng),1);
         for ilp = 1:length(texstrng)
             islegend(ilp) = iscellstr(texstrng{ilp}) && min(strfind([texstrng{ilp}{:}],'input'));
         end
         hleg = texlist(find(islegend,1));
         for tlp = 1:length(cfg.toi)
             hleg.String{tlp} = sprintf('t=% *.*f :%s',fleg,wleg,cfg.toi(tlp),hleg.String{tlp}(end));
         end
     else
         ft_plot_text(0.5, (numel(freqdata{1}.label)+1).*1.2-0.5, sprintf('t=% *.*fs',fleg,wleg,cfg.toi(tlp)), 'interpreter', 'none', 'horizontalalignment', 'right');
     end
     return
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

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'channel'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

k = 1;

% Get physical min/max range of x:
if ischar(cfg.xlim) && strcmp(cfg.xlim,'maxmin')
  xmin = inf;
  xmax = -inf;
  for k = 1:Ndata
    xmin = min(xmin,varargin{k}.(xparam)(1));
    xmax = max(xmax,varargin{k}.(xparam)(end));
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end
cfg.xlim = [xmin xmax];

% Get physical min/max range of y:
if ischar(cfg.ylim) && strcmp(cfg.ylim,'maxmin')
  ymin = inf;
  ymax = -inf;
  for k = 1:Ndata
    ymin = min(ymin,varargin{k}.(yparam)(1));
    ymax = max(ymax,varargin{k}.(yparam)(end));
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end
cfg.ylim = [ymin ymax];

% skip the process for Ndata>1
if nargin > 2, warning('only first data is plotted.');  end
data = varargin{1};

if ~isfield(data, cfg.parameter)
  error('the data does not contain the requested parameter %s', cfg.parameter);
end

% get the selection of the data
tmpcfg           = [];
tmpcfg.latency   = strcmp(xparam,'time')*cfg.xlim + strcmp(yparam,'time')*cfg.ylim;
tmpcfg.frequency = strcmp(xparam,'freq')*cfg.xlim + strcmp(yparam,'freq')*cfg.ylim;
data             = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

dat   = data.(cfg.parameter);
nchan = numel(data.label);


% Get physical min/max range of z:
if ischar(cfg.zlim) && strcmp(cfg.zlim,'maxmin')
  zmin = inf;
  zmax = -inf;
  for k = 1:Ndata
    zmin = min(zmin,min(data.(cfg.parameter)(:)));
    zmax = max(zmax,max(data.(cfg.parameter)(:)));
  end
elseif ischar(cfg.zlim) && strcmp(cfg.zlim,'maxabs')
  zmax = -inf;
  for k = 1:Ndata
    zmax = max(zmax,max(abs(data.(cfg.parameter)(:))));
  end
  zmin = -zmax;
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end
cfg.zlim = [zmin zmax];

% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  set(gcf,'colormap',cfg.colormap);
end

if (isfield(cfg, 'holdfig') && cfg.holdfig==0) || ~isfield(cfg, 'holdfig')
  cla;
  hold on;
end

for k = 1:nchan
  for m = 1:nchan
    if k~=m
      ix  = k;
      iy  = nchan - m + 1;
      % use the convention of the row-channel causing the column-channel
      datamatrix = squeeze(dat(m,k,:,:));
      xvector    = data.(xparam);
      yvector    = data.(yparam);
      
      %-- mask
      if ~isempty(cfg.maskparameter)
          mask = data.(cfg.maskparameter);
          if isfull && cfg.maskalpha == 1
              mask = squeeze(mask(m,k,:,:));
          elseif isfull && cfg.maskalpha ~= 1
              maskl = squeeze(mask(m,k,:,:));
              mask = zeros(size(maskl));
              mask(maskl) = 1;
              mask(~maskl) = cfg.maskalpha;
          end
      end
      
      %-- plot matrix
      if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter)
          nans_mask = ~isnan(datamatrix);
          mask = double(nans_mask);
          ft_plot_matrix(xvector, yvector, datamatrix, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'clim',cfg.zlim, 'box', 'yes','tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
      elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter)
          nans_mask = ~isnan(datamatrix);
          mask = mask .* nans_mask;
          mask = double(mask);
          ft_plot_matrix(xvector, yvector, datamatrix, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'clim',cfg.zlim, 'box', 'yes','tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
      elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter)
          mask = double(mask);
          ft_plot_matrix(xvector, yvector, datamatrix, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'clim',cfg.zlim, 'box', 'yes','tag','cip','highlightstyle',cfg.maskstyle,'highlight', mask)
      else
          ft_plot_matrix(xvector, yvector, datamatrix, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'clim',cfg.zlim, 'box', 'yes','tag','cip')
      end
      
      if k==1,
        % first column, plot scale on y axis
        ft_plot_text( ix.*1.2-0.54,iy.*1.2-0.48,num2str(yvector(1  ),3),'HorizontalAlignment','Right','VerticalAlignment','Middle','Fontsize',cfg.fontsize,'Interpreter','none');
        ft_plot_text( ix.*1.2-0.54,iy.*1.2+0.48,num2str(yvector(end),3),'HorizontalAlignment','Right','VerticalAlignment','Middle','Fontsize',cfg.fontsize,'Interpreter','none');
      end
      if m==nchan,
        % bottom row, plot scale on x axis
        ft_plot_text( ix.*1.2-0.48,iy.*1.2-0.5,num2str(xvector(1  ),3),'HorizontalAlignment','Center','VerticalAlignment','top','Fontsize',cfg.fontsize,'Interpreter','none');
        ft_plot_text( ix.*1.2+0.48,iy.*1.2-0.5,num2str(xvector(end),3),'HorizontalAlignment','Center','VerticalAlignment','top','Fontsize',cfg.fontsize,'Interpreter','none');
      end
    end
  end
end

% add channel labels on grand X and Y axes
maxlablen       = 0;
for k = 1:nchan
  ft_plot_text(0.48,    (nchan + 1 - k).*1.2,     [data.label{k} '  '], 'Interpreter', 'none', 'horizontalalignment', 'right','Fontsize',cfg.fontsize);
  ft_plot_text(k.*1.2,  (nchan + 1)    .*1.2-0.6, ['  ' data.label{k}], 'Interpreter', 'none', 'horizontalalignment', 'left', 'rotation', 90,'Fontsize',cfg.fontsize);
  maxlablen     = max(maxlablen, length(data.label{k}));
end

% add 'from' and 'to' labels
labhight        = sprintf(repmat('\n',1,maxlablen+1));
ft_plot_text(0.25,            (nchan + 1)/1.7, ['\it{from}' labhight], 'rotation', 90,'Fontsize',cfg.fontsize);
ft_plot_text((nchan + 1)/1.7, (nchan + 1)*1.2-0.4, ['\it{to}' labhight],'Fontsize',cfg.fontsize);

axis([-0.2 (nchan+1).*1.2+0.2 0 (nchan+1).*1.2+0.2]);
axis off;

if strcmp(cfg.colorbar,'yes')
  axisSiz   = get(gca, 'Position');
  cylim     = [axisSiz(4) .* (1.*1.2-0.5) ./ ((nchan+1).*1.2+0.2), ...
               axisSiz(4) .* (nchan.*1.2+0.5) ./ ((nchan+1).*1.2+0.2)];
  cylim     = [axisSiz(2) + cylim(1), cylim(2) - cylim(1)];
  cxlim     = [axisSiz(3) .* 1 ./ ((nchan+1).*1.2+0.4), ...
               axisSiz(3) .* ((nchan+0.6).*1.2+0.2) ./ ((nchan+1).*1.2+0.4)]; % slim start at -0.2
  cxlim     = [axisSiz(1)+cxlim(2), cxlim(1)*0.8];
    cxlim   = [cxlim(1), (cxlim(2)*(1+2*(cxlim(2)<0.03)) + 0.03*(3-2*(cxlim(2)<0.03)))/4];   % slightly adjust width to 0.03
  % tag the colorbar so we know which axes are colorbars
  hcb = colorbar('tag', 'ft-colorbar');
%   cbSiz     = get(hcb,'Position');
  set(hcb,'Position', [cxlim(1) cylim(1) cxlim(2) cylim(2)], 'FontSize',cfg.fontsize*0.9);
end

set(gcf, 'color', [1 1 1]);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance

ft_private('remove','silent');
warning(curwarn);

catch ME
    ft_private('remove','silent');
    warning(curwarn);
    rethrow(ME)
end % try