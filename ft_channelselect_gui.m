function label = ft_channelselect_gui(cfg,data)
% 
% FT_CHANNELSELECT_GUI enable to select channels with overlaying topo maps.
% Black dots are the selecting channels.
% 
% Use as
%   [label] = ft_channelselect_gui(cfg, timelock)
% or
%   [label] = ft_channelselect_gui(cfg, freq)
% or
%   [label] = ft_channelselect_gui(cfg)
% 
% 'cfg' options are basically same as 'ft_topoplotER'.
% 
% left click:   toggle
% right click:  add    (ctrl+click)
% middle click: remove (shift+click)
% 
% press 'q' or 'esc' to quit
%  (press 'q' to close the figure)
% 
% Example 1: 
%     tmpcfg = [];
%     tmpcfg.latency = [0.160 0.180];
%     tmpcfg.layout = 'neuromag306mag.lay';
%     tmpcfg.colormap = 'jet';
%     tmpcfg.selectchannel = {'MEG243*','MEG244*','MEG241*','MEG232*','MEG251*','MEG252*'};
%     cfg.channel = ft_channelselect_gui(tmpcfg,avg);
% 
% Example 2: 
%     tmpcfg = [];
%     tmpcfg.layout = 'neuromag306cmb.lay';
%     cfg.channel = ft_channelselect_gui(tmpcfg);
% 
% Output channels are extended to all channels at the same positions based
% on data.grad.
% If you want to select only magnetmeter/gradiometer, reject the other
% channels.
% 
% Example:
%     %-- select megnetmeter
%     cfg.channel = [cfg.channel; {'-MEGGRAD'}];
%     %-- select gradiometer
%     cfg.channel = [cfg.channel; {'-MEGMAG'}];
% 
% If layout of combined channels is indicated and input data is not
% combined, ft_combineplanar is automatically applied.

% using: fieldtrip, ft_private(ft_platform_supports), cellstrfind,
%        findarray

% 20160822 Yuasa
% 20160906 Yuasa: rename 'latency' -> 'xlim', 'selectchannel' -> 'highlightchannel'
%                 automatically apply ft_combineplanner
% 20170517 Yuasa: enable without data
% 20180404 Yuasa: enable TFR data
% 20181225 Yuasa: bug fix
% future -> add culculate button,  selectchannel on combined data

ft_defaults

if isfield(cfg,'latency') && ~isfield(cfg,'xlim')
    cfg.xlim  = cfg.latency;
    warning('set cfg.xlim based on cfg.latency');
end
if isfield(cfg,'selectchannel') && ~isfield(cfg,'highlightchannel')
    cfg.highlightchannel  = cfg.selectchannel;
    warning('set cfg.highlightchannel based on cfg.selectchannel');
end

cfg.colorbar          = ft_getopt(cfg, 'colorbar',      'no');
cfg.colormap          = ft_getopt(cfg, 'colormap',      'jet');
cfg.shading           = ft_getopt(cfg, 'shading',       'flat');
cfg.comment           = ft_getopt(cfg, 'comment',       'no');
cfg.fontsize          = ft_getopt(cfg, 'fontsize',      11);
cfg.baseline          = ft_getopt(cfg, 'baseline',      'no'); %to avoid warning in timelock/freqbaseline
cfg.interactive       = ft_getopt(cfg, 'interactive',   'no');
cfg.marker            = ft_getopt(cfg, 'marker',        'on');
cfg.markersymbol      = ft_getopt(cfg, 'markersymbol',  'o');
cfg.markercolor       = ft_getopt(cfg, 'markercolor',   [0 0 0]);
cfg.markersize        = ft_getopt(cfg, 'markersize',    5);
cfg.markerfontsize    = ft_getopt(cfg, 'markerfontsize', 11);
cfg.highlight         = ft_getopt(cfg, 'highlight',     'on');
cfg.highlightchannel  = ft_getopt(cfg, 'highlightchannel',  '', 1); % highlight may be 'on', making highlightchannel {} meaningful
cfg.highlightsymbol   = ft_getopt(cfg, 'highlightsymbol',   '.');
cfg.highlightcolor    = ft_getopt(cfg, 'highlightcolor',    [0 0 0]);
cfg.highlightsize     = ft_getopt(cfg, 'highlightsize',     20);
cfg.highlightfontsize = ft_getopt(cfg, 'highlightfontsize', 11);
cfg.labeloffset       = ft_getopt(cfg, 'labeloffset',       0.005);
cfg.channel           = ft_getopt(cfg, 'channel',           'all');
cfg.figurename        = ft_getopt(cfg, 'figurename',        []);
cfg.interpolatenan    = ft_getopt(cfg, 'interpolatenan',    'yes');
cfg.trials            = ft_getopt(cfg, 'trials',        'all', 1);

hasdata     = nargin > 1;
haslay      = isfield(cfg,'layout');
%-- make sample data
if ~hasdata
    assert(haslay, 'please specify layout or data');
    layout      = ft_prepare_layout(cfg);
else
    layout      = ft_prepare_layout(cfg,data);
end
layout.label(match_str(layout.label,{'COMNT';'SCALE'})) = [];
if ~hasdata
    cfg = ft_checkconfig(cfg, 'renamed',    {'zparam', 'parameter'});
    cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
    
    data        = [];
    data.time   = 0;
    data.label  = layout.label;
    data.dimord = 'chan_time';
    data.(cfg.parameter) = zeros(length(data.label),length(data.time));
end
hasgrad       = isfield(data,'grad');

%-- combine planar
if hasgrad && isfield(data.grad,'type')
    iscombined = ~isempty(regexp(data.grad.type,'.+_combined$', 'once'));   % grad.type = '..._combined'
else
    iscombined = ~isempty(cellstrfind(data.label,'*+*'));           % data.label = {'...+...'}
end
isoutputcmb   = ~isempty(cellstrfind(layout.label,'*+*'));          % layout.label = {'...+...'}
    
if isoutputcmb && ~iscombined
    if hasgrad,  orig_grad   = data.grad;    end
    tmpcfg              = [];
    tmpcfg.method       = 'sum';
    tmpcfg.demean       = 'no';
    data        = ft_combineplanar(tmpcfg, data);
    if hasgrad
        if ~haslay
            data.grad       = reducegrad(data.label, data.grad);      % pick up combined channel
            cfg.layout      = ft_prepare_layout(cfg,data);
        end
        %-- restore grad
        data.combinegrad    = data.grad;
        data.grad           = orig_grad;
    end
end

%-- prepare layout
if ~isfield(cfg,'layout')
    cfg.layout      = ft_prepare_layout(cfg,data);
end

%-- draw topography
figure;
cfgtp=ft_topoplotER(cfg,data);
htp = gcf;

%-- For Highlight (channel-selection)
    highlightchansel = []; % used for remembering selection of channels
    for icell = 1:length(cfgtp.highlight)
      if ~strcmp(cfgtp.highlight{icell},'off')
        highlightchansel   = [highlightchansel; match_str(data.label,ft_channelselection(cfgtp.highlightchannel{icell}, cfgtp.layout.label))];
      end
    end
    cfg.highlightchannel = data.label(highlightchansel);
    cfg.data             = data;
 
  %-- add the channel position information to the figure
  %-- this section is important for select_channel_single/multiple
  info        = guidata(htp);
    labellist = true(length(cfgtp.layout.label),1);       % reject 'COMNT''SCALE'
    labellist(cellstrfind(cfgtp.layout.label,{'COMNT','SCALE'})) = false;
  info.x      = cfgtp.layout.pos(labellist,1);
  info.y      = cfgtp.layout.pos(labellist,2);
  info.label  = cfgtp.layout.label(labellist);
  info.cfg    = cfg;
  info.delete = 0;
  guidata(htp, info);
  
  %-- multiple selection
  set(htp, 'KeyPressFcn',           @keyboard_cb);
  set(htp, 'WindowButtonDownFcn',   @startselectrange);
  set(htp, 'WindowButtonUpFcn',     @decidechannel);
  set(htp, 'WindowButtonMotionFcn', @selectrange);

  
  %-- output
  if nargout > 0
    uiwait;
    
    if ishandle(htp)
      info  = guidata(htp);
      label = info.cfg.highlightchannel;
      
      %-- extend select channel to all mag & planer
      if hasgrad
          cur_label = label;
          label    = [];
          for selich = 1:length(cur_label)
              if isfield(data,'combinegrad')
                  %-- find channels at the same position to combined sensor
                  gradidx = cellstrfind(data.combinegrad.label,cur_label{selich});
                  chidx = findarray(data.grad.chanpos,data.combinegrad.chanpos(gradidx,:));
              else
                  %-- find planars at the same position to magnetometer
                  gradidx = cellstrfind(data.grad.label,cur_label{selich});
                  chidx = findarray(data.grad.chanpos,data.grad.chanpos(gradidx,:));
              end
              label = [label; data.grad.label(chidx)];
          end
      end
      label = unique(label);
      
      if info.delete
          delete(htp);
      end
    else
      label = [];
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of a single channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_channel_single(~)
htp = gcf;
%-- get last-used-mouse-button
lastmousebttn = get(gcf,'selectiontype');

pos = get(gca, 'CurrentPoint');
pos = pos(1,1:2);

info  = guidata(gcf);
x     = info.x;
y     = info.y;
label = info.label;
cfg    = info.cfg;

%-- compute a tolerance measure
distance = sqrt(abs(sum([x y]'.*[x y]',1)));
distance = triu(distance, 1);
distance = distance(:);
distance = distance(distance>0);
distance = median(distance);
tolerance = 0.3*distance;

%-- compute the distance between the clicked point and all channels
dx = x - pos(1);
dy = y - pos(2);
dd = sqrt(dx.^2 + dy.^2);
[d, i] = min(dd);
if d<tolerance
  label = label{i};
  fprintf('channel "%s" selected\n', label);
  
  toggle_ch = match_str(cfg.highlightchannel, label);
  if isempty(toggle_ch) && ~strcmp(lastmousebttn,'extend')
      cfg.highlightchannel = cat(1, cfg.highlightchannel, label);
  elseif ~isempty(toggle_ch) && ~strcmp(lastmousebttn,'alt')
      cfg.highlightchannel(toggle_ch) = [];
  end
  
  %-- clear & replot
  clf(htp); figure(htp);
  ft_topoplotER(cfg,cfg.data);
  info.cfg = cfg;
  guidata(htp, info);
else
  label = {};
  fprintf('no channel selected\n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of multiple channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_channel_multiple(range) % last input is context menu label, see ft_select_range
htp = gcf;
%-- get last-used-mouse-button
lastmousebttn = get(gcf,'selectiontype');
  
info   = guidata(htp);
x      = info.x(:);
y      = info.y(:);
label  = info.label(:);
cfg    = info.cfg;

%-- determine which channels ly in the selected range
select = false(size(label));
for i=1:size(range,1)
  select = select | (x>=range(i, 1) & x<=range(i, 2) & y>=range(i, 3) & y<=range(i, 4));
end
label  = label(select);
if isempty(label)
  fprintf('no channel selected\n');
else
  fprintf('channel %s selected\n', sprintf('"%s"',label{:}));
  
  switch lastmousebttn
      case 'normal' % left click: toggle
          [toggle_ch1, toggle_ch2] = match_str(cfg.highlightchannel, label);
          cfg.highlightchannel(toggle_ch1) = [];
          label(toggle_ch2) = [];
          cfg.highlightchannel = cat(1, cfg.highlightchannel, label);
      case 'alt' % right click: add
          [toggle_ch1, toggle_ch2] = match_str(cfg.highlightchannel, label);
          label(toggle_ch2) = [];
          cfg.highlightchannel = cat(1, cfg.highlightchannel, label);
      case 'extend' % shift-click: remove
          [toggle_ch1, toggle_ch2] = match_str(cfg.highlightchannel, label);
          cfg.highlightchannel(toggle_ch1) = [];
  end
  %-- clear & replot
  clf(htp); figure(htp);
  ft_topoplotER(cfg,cfg.data);
  info.cfg = cfg;
  guidata(htp, info);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of multiple channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startselectrange(handle,~) % last input is context menu label, see ft_select_range
[handle, p, UserData, selecting] = prepareselect(handle);

%-- get last-used-mouse-button
lastmousebttn = get(gcf,'selectiontype');

    switch lastmousebttn
      case {'normal', 'alt', 'extend'} % left click, right click, shift-click
          % ~multiple
            %-- start with a new selection
            delete(UserData.box(ishandle(UserData.box)));
            UserData.range = [];
            UserData.box   = [];
          
          %-- add a new selection range
          UserData.range(end+1,1:4) = nan;
          UserData.range(end,1) = p(1);
          UserData.range(end,3) = p(2);
          
          %-- add a new selection box
          xData = [nan nan nan nan nan];
          yData = [nan nan nan nan nan];
          UserData.box(end+1) = line(xData, yData);
    end

%-- put the modified selections back into the figure
if ishandle(handle)
  setappdata(handle, 'select_range_m', UserData);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of multiple channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectrange(handle,~) % last input is context menu label, see ft_select_range
[handle, p, UserData, selecting] = prepareselect(handle);

if selecting
    %-- update the selection box
        x1 = UserData.range(end,1);
        x2 = p(1);
        y1 = UserData.range(end,3);
        y2 = p(2);
    
    xData = [x1 x2 x2 x1 x1];
    yData = [y1 y1 y2 y2 y1];
    set(UserData.box(end), 'xData', xData);
    set(UserData.box(end), 'yData', yData);
    set(UserData.box(end), 'Color', [0 0 0]);
    %set(UserData.box(end), 'EraseMode', 'xor');
    set(UserData.box(end), 'LineStyle', '--');
    set(UserData.box(end), 'LineWidth', 1.5);
    set(UserData.box(end), 'Visible', 'on');
    
else
        set(handle, 'Pointer', 'crosshair');
end

%-- put the modified selections back into the figure
if ishandle(handle)
  setappdata(handle, 'select_range_m', UserData);
end

    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of multiple channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decidechannel(handle,~) % last input is context menu label, see ft_select_range
[handle, p, UserData, selecting] = prepareselect(handle);

%-- get last-used-mouse-button
lastmousebttn = get(gcf,'selectiontype');

    switch lastmousebttn
      case {'normal', 'alt', 'extend'} % left click, right click, shift-click

        if selecting
          %-- select the other corner of the box
          UserData.range(end,2) = p(1);
          UserData.range(end,4) = p(2);
        end
                
        if ~isempty(UserData.range)
          %-- ensure that the selection is sane
          if diff(UserData.range(end,1:2))<0
            UserData.range(end,1:2) = UserData.range(end,[2 1]);
          end
          if diff(UserData.range(end,3:4))<0
            UserData.range(end,3:4) = UserData.range(end,[4 3]);
          end
        end

          %-- finish -> select channel
          if diff(UserData.range(end,1:2)) + diff(UserData.range(end,1:2)) < 10^-3
              callback = 'select_channel_single';
          else
              callback = 'select_channel_multiple';
          end
          feval(callback, UserData.range);
          delete(UserData.box(ishandle(UserData.box)));
          UserData.range = [];
          UserData.box   = [];
          set(handle, 'Pointer', 'crosshair');
        
    end

%-- put the modified selections back into the figure
if ishandle(handle)
  setappdata(handle, 'select_range_m', UserData);
end
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of multiple channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [handle, p, UserData, selecting] = prepareselect(handle)
ft_private('silent');
p = handle;
if ft_platform_supports('graphics_objects')
 while ~isa(handle, 'matlab.ui.Figure')
    handle = p;
    p = get(handle, 'parent');
 end
else
    while ~isequal(p, 0) 
      handle = p;
      p = get(handle, 'parent');
    end
end
ft_private('remove','silent');

if ishandle(handle)
  UserData = getappdata(handle, 'select_range_m');
else
  UserData = [];
end

if isempty(UserData)
  UserData.range = []; % this is a Nx4 matrix with the selection range
  UserData.box   = []; % this is a Nx1 vector with the line handle
end

p = get(gca, 'CurrentPoint');
p = p(1,1:2);

abc = axis;
xLim = abc(1:2);
yLim = abc(3:4);

%-- limit cursor coordinates
if p(1)<xLim(1), p(1)=xLim(1); end;
if p(1)>xLim(2), p(1)=xLim(2); end;
if p(2)<yLim(1), p(2)=yLim(1); end;
if p(2)>yLim(2), p(2)=yLim(2); end;

%-- determine whether the user is currently making a selection
selecting = numel(UserData.range)>0 && any(isnan(UserData.range(end,:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboard_cb(handle, eventdata)
ft_private('silent');
if (isempty(eventdata) && ft_platform_supports('matlabversion',-Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  %-- determine the key that corresponds to the uicontrol element that was activated
  key = get(handle, 'userdata');
else
  %-- determine the key that was pressed on the keyboard
  key = eventdata.Key;
end
ft_private('remove','silent');

htp = gcf;

switch key
    case 'q'
        info        = guidata(htp);
        info.delete = 1;
        guidata(htp,info);
        uiresume;
    case 'escape'
        uiresume;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = reducegrad(label, grad)
labelidx = match_str(grad.label, label);
rmidx    = 1:length(grad.label);        rmidx(labelidx) = [];

if isfield(grad,'chanori'),  grad.chanori(rmidx,:) = []; end;
if isfield(grad,'chanpos'),  grad.chanpos(rmidx,:) = []; end;
if isfield(grad,'chantype'), grad.chantype(rmidx,:) = []; end;
if isfield(grad,'chanunit'), grad.chanunit(rmidx,:) = []; end;
if isfield(grad,'newtype'),  grad.newtype(rmidx,:) = []; end;
if isfield(grad,'newunit'),  grad.newunit(rmidx,:) = []; end;
grad.label(rmidx,:) = [];

