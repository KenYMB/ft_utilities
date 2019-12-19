function ft_plot_grangernet(cfg,data)
% FT_PLOT_GRANGERNET visualize the Granger causality network.
% The line of granger network is judged by the max peak causality along
% with frequencies in each time point.
% 
% Syntax:
%   ft_plot_grangernet(cfg.data)
% 
% The cfg can have the following options:
%   cfg.parameter   = string, the functional parameter to be plotted (default = 'grangerspctrm')
%   cfg.toi         = number, time of interest
%   cfg.foilim      = [begin end], frequency band of interest
%                     'all' or [fmin fmax] (default = 'all')
%   cfg.channel     = list of channels to be included for the plotting (default = 'all')
%   cfg.location    = 2 x channels vector, Location of the sites
%                     'auto' or [xpos; ypos] (default = 'auto')
%   cfg.threshold   = number, threshold for the plotting (default = '10%')
%   (cfg.fontsize    = font size of this figure (default = 9))
% 
% 
% See also: BS_PLOT_GRANGERNET, FT_CONNECTIVITYPLOTTFR, BS2FT_GRANGERANALYSIS, BS_GRANGER.

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 14-Sep-2007 10:38:30$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% Using: draw_layout (, make_layout, arrow, textoval, ellipse),
%        nearlyeq, nearvalue

% 20170329 Yuasa: modified "ganetwork" to use fieldtrip format

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

%-- do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance varargin
ft_preamble trackconfig

%-- the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'foi', 'foilim'});

%-- parameter settings% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel',      'all');
cfg.parameter  = ft_getopt(cfg, 'parameter',    'grangerspctrm');
cfg.toi        = ft_getopt(cfg, 'toi',      []);
cfg.foilim     = ft_getopt(cfg, 'foilim',   'all');
cfg.location   = ft_getopt(cfg, 'location', 'auto');
cfg.threshold  = ft_getopt(cfg, 'threshold','auto');
    if strcmp(cfg.threshold,'auto'),    cfg.threshold = '10%';  end
cfg.fontsize   = ft_getopt(cfg, 'fontsize',     9);

%-- set private functions
ft_private('set','silent')

try
    
  %-- check if the input data is valid for this function
  data    = ft_checkdata(data, 'datatype', {'timelock', 'freq'});
  dimord      = getdimord(data, cfg.parameter);
  dimtok      = tokenize(dimord, '_');

  %-- convert into the the supported dimord
  dim_chan           = find(strcmp(dimtok,'chan'));
  dim_chancmb        = find(strcmp(dimtok,'chancmb'));
  dim_freq           = find(strcmp(dimtok,'freq'));
  dim_time           = find(strcmp(dimtok,'time'));
    
  if length(dim_chan)==2 && length(dim_freq)==1 && length(dim_time)==1
    %-- that's ok
    data.(cfg.parameter)   = permute(data.(cfg.parameter), [dim_chan dim_chancmb dim_freq dim_time]);
  elseif length(dim_chancmb)==1 && length(dim_freq)==1 && length(dim_time)==1
    %-- convert into 'chan_chan_freq_time'
    data.(cfg.parameter)   = permute(data.(cfg.parameter), [dim_chan dim_chancmb dim_freq dim_time]);
    data = ft_checkdata(data, 'cmbrepresentation', 'full');
  else
    error('the data should have a dimord of %s or %s', 'chan_chan_freq_time', 'chancmb_freq_time');
  end
  
  % set x_y param
  Nchans = length(data.label);
  gcpowmax  = round(nanmax(data.(cfg.parameter)(:)),3,'significant');
  gcpowmin  = 0;
  
  %-- parameter settings
  if ischar(cfg.channel) && strcmp(cfg.channel,'all')
        cfg.channel   = 1:Nchans;
  else
        cfg.channel   = round(cfg.channel);
        cfg.channel(cfg.channel<1) = [];
        cfg.channel(cfg.channel>Nchans) = [];
  end
  if isempty(cfg.toi),  cfg.toi = data.time(1);
  else                  cfg.toi = nearvalue(data.time,cfg.toi(1));
  end
  if ischar(cfg.foilim) && strcmp(cfg.foilim,'all')
        cfg.foilim   = data.freq([1, end]);
  else  cfg.foilim   = nearvalue(data.freq,cfg.foilim(1:2));
  end
  if ischar(cfg.location) && strcmp(cfg.location,'auto')
        cfg.location   = prepare_layout(Nchans);
  elseif size(cfg.location,2) < Nchans
        cfg.location   = prepare_layout(Nchans);
        warning('location data is lacked. use default layout.');
  else  cfg.location   = cfg.location(1:2,1:Nchans);
  end
  if ischar(cfg.threshold) && strcmp(cfg.threshold(end),'%')
        cfg.threshold    = round(gcpowmax * str2double(cfg.threshold(1:end-1)) / 100, 3, 'significant');
  else  cfg.threshold    = max(min(cfg.threshold,gcpowmax),gcpowmin);
  end
  
  
  %-- main
gch = figure('Name','Granger Causality Network','NumberTitle','off','Menubar','none','ToolBar','none');
  gpa = plot_gc_layout(cfg,data);
  gpa.Tag = 'Main';
  
  setappdata(gch,'cfg',cfg);
  setappdata(gch,'data',data);
  
  axpos = gpa.Position;
  axrsp = 1-(axpos(1)+axpos(3));
  axusp = 1-(axpos(2)+axpos(4));
  htl = uicontrol('Style','text','Tag','Tlabel',...
                  'Units','normalized','Position',[axpos(1)-0.09 axpos(2)/3 0.08 axpos(2)/3],...
                  'String','Time','HorizontalAlignment','right');
  hsl = uicontrol('Style','slider','Tag','Tslider',...
                  'Units','normalized','Position',[axpos(1) axpos(2)/3 axpos(3) axpos(2)/3],...
                  'Min',data.time(1),'Max',data.time(end),'SliderStep',[1/(length(data.time)-1) 2/(length(data.time)-1)],...
                  'Value',cfg.toi, 'Callback',@timesld);
  hfl = uicontrol('Style','text','Tag','Flabel',...
                  'Units','normalized','Position',[(axpos(1)+axpos(3))+axrsp/8 (axpos(2)+axpos(4))+axusp/20 axrsp*3/4 axusp/2],...
                  'String','Freq');
  hfd = uicontrol('Style','slider','Tag','FMINslider',...
                  'Units','normalized','Position',[(axpos(1)+axpos(3))+axrsp/6 axpos(2) axrsp/3 axpos(4)],...
                  'Min',data.freq(1),'Max',data.freq(end),'SliderStep',[1/(length(data.freq)-1) 5/(length(data.freq)-1)],...
                  'Value',cfg.foilim(1), 'Callback',@freqsld);
  hfu = uicontrol('Style','slider','Tag','FMAXslider',...
                  'Units','normalized','Position',[(axpos(1)+axpos(3))+axrsp/2 axpos(2) axrsp/3 axpos(4)],...
                  'Min',data.freq(1),'Max',data.freq(end),'SliderStep',[1/(length(data.freq)-1) 5/(length(data.freq)-1)],...
                  'Value',cfg.foilim(2), 'Callback',@freqsld);
  hthr = uicontrol('Style','pushbutton','Tag','THREb',...
                  'Units','normalized','Position',[axpos(1)/5 (axpos(2)+axpos(4))+axusp/10 axpos(1)*4/5 axusp*3/5],...
                  'String','Threshold','Callback',@threshedit);
  hthsl = uicontrol('Style','slider','Tag','THREslider',...
                  'Units','normalized','Position',[axpos(1)/2 axpos(2) axrsp/3 axpos(4)],...
                  'Min',gcpowmin,'Max',gcpowmax,'SliderStep',[round(gcpowmax/20,1,'significant') round(gcpowmax/5,1,'significant')],...
                  'Value',cfg.threshold, 'Callback',@threshsld);
  
  
  %-- do the general cleanup and bookkeeping at the end of the function
  ft_postamble debug
  ft_postamble trackconfig
  ft_postamble previous varargin
  ft_postamble provenance
  
  ft_private('remove','silent');
  
catch ME
  %-- do the general cleanup and bookkeeping at the end of the function
  ft_postamble debug
  ft_postamble trackconfig
  ft_postamble previous varargin
  ft_postamble provenance
  
  ft_private('remove','silent');
  
  rethrow(ME);
end

end%bs_plot_grangernet

function hax = plot_gc_layout(cfg,data)

%-- main
[dag, ispeak] = gc_data_processing(cfg,data);

%-- plot figure
tstr = sprintf('Time: %.2g s, Frequency: %.0f - %.0f Hz, Threshold: %.2g',cfg.toi,cfg.foilim(1),cfg.foilim(2),cfg.threshold);
if ~ispeak
    tstr = sprintf('No peak found. %s',tstr);
end
%-- set parameters for draw
Nchans = length(data.label);
circle = zeros(Nchans,1);
circle(cfg.channel) = 1; 
left   = cfg.location(1,:);
right  = cfg.location(2,:);
%-- plot main
cla;
draw_layout(dag,data.label,circle,left,right,'gc');
hax = gca;
title(tstr,'Tag','Title');

set(hax.Children.findobj('Type','text'),'ButtonDownFcn', @chansel);
set(hax.Children.findobj('Type','patch','Tag',''),'ButtonDownFcn', @chansel);

end%plot_gc_layout

function [dag, ispeak] = gc_data_processing(cfg,data)

Nchans = length(data.label);
toiidx = nearlyeq(data.time, cfg.toi);
foiidx = nearlyeq(data.freq,cfg.foilim(1)):nearlyeq(data.freq,cfg.foilim(2));

ispeak  = true;
dag     = zeros(Nchans,Nchans);
dag2    = dag;
 
for ich = 1:Nchans
  for jch = [1:(ich-1) (ich+1):Nchans]
    peak=findpeaks(squeeze(data.(cfg.parameter)(ich,jch,foiidx,toiidx)));       % need SignalProcessingToolbox
    if ~isempty(peak) && max(peak)>cfg.threshold
        dag(ich,jch)=1;
    end
  end
end

%-- set channels of interests
chanary = false(1,Nchans);
chanary(cfg.channel) = true;
dag(~chanary,:) = 0;
dag(:,~chanary) = 0;

if max(dag(:)==dag2(:))
%-- if no peak is found
    ispeak  = false;
    for ich = 1:Nchans
      for jch = [1:(ich-1) (ich+1):Nchans]
        powmax = max(squeeze(data.(cfg.parameter)(ich,jch,foiidx,toiidx)));
        if powmax>cfg.threshold
            dag(ich,jch)=1;
        end
      end
    end
    %-- remove channels of no-interests
    dag(~chanary,:) = 0;
    dag(:,~chanary) = 0;
end

end%gc_data_processing

function layout_data = prepare_layout(num)
% make default layout map

%-- set square or rectangle
nx = ceil(sqrt(num));
lstclmn  = mod(num,nx);
lstclmn2 = mod(num,nx+1);
if lstclmn && (~lstclmn2 || lstclmn2 >= lstclmn) 
  nx = nx + 1;
end
ny = ceil(num/nx);

posx = 0.9/(nx+1)*(1:nx)+0.05;
posy = 0.9/(ny+1)*sort([1:2:ny (2:2:ny)+0.1])+0.05;

[meshx, meshy] = meshgrid(posx, posy);
if nx >= 3
  meshy(:,1:2:nx) = meshy(:,1:2:nx) - (meshy(1,1)-0.05)*0.5/2;
  meshy(:,2:2:nx) = meshy(:,2:2:nx) + (meshy(1,1)-0.05)*0.5/2;
end
if ny >= 3
  meshx(1:2:ny,:) = meshx(1:2:ny,:) - (meshx(1,1)-0.05)*0.2/2;
  meshx(2:2:ny,:) = meshx(2:2:ny,:) + (meshx(1,1)-0.05)*0.2/2;
end

layout_data   = [reshape(meshx,1,[]); reshape(meshy,1,[])];

end%prepare_layout

function c = num2cellstr(a)
    c = cell(size(a));
    for i=1:numel(a)
        c{i} = num2str(a(i));
    end 
end%num2cellstr


%----------------- GUI sub function -----------------------%
function timesld(hsl, cld)
    gch = hsl.Parent;
    
    cfg      = getappdata(gch,'cfg');
    data     = getappdata(gch,'data');
    
    toiidx      = nearlyeq(data.time,hsl.Value);
    hsl.Value   = data.time(toiidx);
    if cfg.toi ~= hsl.Value
        cfg.toi = hsl.Value;
        setappdata(gch,'cfg',cfg);
        
        gpa = plot_gc_layout(cfg,data);
    end

end%timesld

function freqsld(hfreq, cld)
    gch = hfreq.Parent;
    hfu = gch.Children.findobj('Tag','FMAXslider');
    hfd = gch.Children.findobj('Tag','FMINslider');
    
    cfg      = getappdata(gch,'cfg');
    data     = getappdata(gch,'data');
        
    foiidx      = nearlyeq(data.freq,hfreq.Value);
    switch hfreq.Tag
        case 'FMAXslider'
            hfreq.Value = max(hfd.Value+2,data.freq(foiidx));
        case 'FMINslider'
            hfreq.Value = min(hfu.Value-2,data.freq(foiidx));
    end
    if max(cfg.foilim ~= [hfd.Value hfu.Value])
        cfg.foilim = [hfd.Value hfu.Value];
        setappdata(gch,'cfg',cfg);
        
        plot_gc_layout(cfg,data);
    end

end%freqsld

function threshedit(hthr, cld)
    gch = hthr.Parent;
    hthsl = gch.Children.findobj('Tag','THREslider');
    
    cfg      = getappdata(gch,'cfg');
    data     = getappdata(gch,'data'); 
    
    hthr.Value = str2double(cell2mat(inputdlg('Threshold?','',1,num2cellstr(cfg.threshold))));
    if isnan(hthr.Value), hthr.Value = cfg.threshold;  end
    if hthr.Value < hthsl.Min, hthr.Value = hthsl.Min;  end
    if hthr.Value > hthsl.Max, hthr.Value = hthsl.Max;  end

    if cfg.threshold ~= hthr.Value
        cfg.threshold = hthr.Value;
        hthsl.Value = cfg.threshold;
        setappdata(gch,'cfg',cfg);
        
        plot_gc_layout(cfg,data);
    end

end%threshedit

function threshsld(hthsl, cld)
    gch = hthsl.Parent;
    
    cfg      = getappdata(gch,'cfg');
    data     = getappdata(gch,'data');
    
    hthsl.Value = round(hthsl.Value,3,'significant');

    if cfg.threshold ~= hthsl.Value
        cfg.threshold = hthsl.Value;
        setappdata(gch,'cfg',cfg);
        
        plot_gc_layout(cfg,data);
    end

end%threshsld

function chansel(hobj, cld)
    gpa = hobj.Parent;
    gch = gpa.Parent;
    Mpos = gpa.CurrentPoint;
    
    cfg      = getappdata(gch,'cfg');
    data     = getappdata(gch,'data');
    
    Mdist    = mean((cfg.location - repmat(Mpos(1,1:2)',1,size(cfg.location,2))).^2,1);
    selchidx = find(Mdist == min(Mdist),1);
    
    if max(cfg.channel == selchidx)
        cfg.channel(cfg.channel == selchidx) = [];
    else
        cfg.channel = sort([cfg.channel selchidx]);
    end
    
    setappdata(gch,'cfg',cfg);
    plot_gc_layout(cfg,data);

end%chansel
% [EOF]