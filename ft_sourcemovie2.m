function [cfg, M] = ft_sourcemovie2(cfg, source, source2)

% FT_SOURCEMOVIE displays the source reconstruction on a cortical mesh
% and allows the user to scroll through time with a movie
% This function is updated from original
%
% Use as
%   ft_sourcemovie(cfg, source)
% where the input source data is obtained from FT_SOURCEANALYSIS and cfg is
% a configuratioun structure that should contain
%
%  cfg.funparameter    = string, functional parameter that is color coded (default = 'pow')
%  cfg.maskparameter   = string, functional parameter that is used for opacity (default = [])
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEPLOT, FT_SOURCEINTERPOLATE
% 
% Usage example
%     cfg = [];
%     cfg.projectmom  = 'yes';
%     sdMEG           = ft_sourcedescriptives(cfg, sourceMEG);
% 
%     cfg = [];
%     cfg.funparameter  = 'pow';
%     cfg.maskparameter = 'pow';
%     cfg.xlim          = [0 0.5];
%     cfg.zlim          = [0 50];
%     cfg.alim          = [10 60];
%     cfg.colormap      = 'hot';
%     ft_sourcemovie(cfg,sdMEG);

% Undocumented options:
%   cfg.parcellation
%   cfg.atlas
%   cfg.queryrange

% Copyright (C) 2011-2015, Robert Oostenveld
% Copyright (C) 2012-2014, Jorn Horschig
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

% 20160912 Yuasa: fix for MATLAB R2016a
%                 add right-click & drag to rotate
% 20160913 Yuasa: bug fix
% 20170321 Yuasa: minor fix
% 20170713 Yuasa: enable to modity alim
% 20170809 Yuasa: enable cfg.altas option
% 20170809 Yuasa: rename to 'ft_sourcemovie2'

% Using: ft_private(ft_platform_supports, parameterselection, getdimord, [getdimsiz], intersect_line, [normals], atlas_lookup)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar source
ft_preamble provenance source
ft_preamble trackconfig

ft_private('silent');

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input argument or can be read from disk
hassource2 = exist('source2', 'var');

% check if the input data is valid for this function
if isfield(cfg, 'atlas') && ~isempty(cfg.atlas)
  % the atlas lookup requires the specification of the coordsys
    source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes', 'hascoordsys', 'yes');
else
    source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',	 {'zparam',    'cfg.funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'parameter', 'cfg.funparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'mask',      'maskparameter'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'olim',      'alim'});
cfg = ft_checkconfig(cfg, 'renamed',	 {'clim',      'zlim'});

% these are not needed any more, once the source structure has a proper dimord
% cfg = ft_checkconfig(cfg, 'deprecated', 'xparam');
% cfg = ft_checkconfig(cfg, 'deprecated', 'yparam');

% get the options
xlim              = ft_getopt(cfg, 'xlim');
ylim              = ft_getopt(cfg, 'ylim');
zlim              = ft_getopt(cfg, 'zlim');
olim              = ft_getopt(cfg, 'alim');                           % don't use alim as variable name
cfg.xparam        = ft_getopt(cfg, 'xparam');                         % default is dealt with below
cfg.yparam        = ft_getopt(cfg, 'yparam');                         % default is dealt with below
cfg.funparameter  = ft_getopt(cfg, 'funparameter');
cfg.maskparameter = ft_getopt(cfg, 'maskparameter');
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');
cfg.title         = ft_getopt(cfg, 'title',         '');
cfg.parcellation  = ft_getopt(cfg, 'parcellation');
cfg.atlas         = ft_getopt(cfg, 'atlas');
cfg.queryrange    = ft_getopt(cfg, 'queryrange',          3);
cfg.colormap      = ft_getopt(cfg, 'colormap','default');
if ~isfield(source, 'tri')
    cfg.surfdownsample = ft_getopt(cfg, 'surfdownsample', 1);
    cfg.surffile       = ft_getopt(cfg, 'surffile', 'surface_white_both.mat'); % use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
    cfg.surftrans      = ft_getopt(cfg, 'surftrans');
    cfg.surfinflated   = ft_getopt(cfg, 'surfinflated');
    cfg.sphereradius   = ft_getopt(cfg, 'sphereradius',  []);
    cfg.projvec        = ft_getopt(cfg, 'projvec',       1);
    cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
    cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); % or max
    cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);
    cfg.projmethod     = ft_getopt(cfg, 'projmethod',    'nearest');
    cfg.distmat        = ft_getopt(cfg, 'distmat',       []);
end

% select the functional and the mask parameter
cfg.funparameter  = parameterselection(cfg.funparameter, source);
cfg.maskparameter = parameterselection(cfg.maskparameter, source);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1};  end
try, cfg.maskparameter = cfg.maskparameter{1}; end

dimord = getdimord(source, cfg.funparameter);
dimtok = tokenize(dimord, '_');

if isempty(cfg.xparam) && numel(dimtok)>1
  cfg.xparam = dimtok{2};
end

if isempty(cfg.yparam) && numel(dimtok)>2
  cfg.yparam = dimtok{3};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~hassource2
  fun = getsubfield(source, cfg.funparameter);
elseif hassource2 && isfield(source2, 'pos'),
  fun  = getsubfield(source, cfg.funparameter);
  fun2 = getsubfield(source2, cfg.funparameter);
elseif hassource2
  % assume the first data argument to be a parcellation, and the second a parcellated structure
  tmp = getsubfield(source2, cfg.funparameter);
  siz = [size(tmp) 1];
  fun = zeros([size(source.pos, 1), siz(2:end)]);
  parcels      = source.(cfg.parcellation);
  parcelslabel = source.([cfg.parcellation,'label']);
  for k = 1:numel(source2.label)
    sel = match_str(source.([cfg.parcellation,'label']), source2.label{k});
    if ~isempty(sel)
      sel = source.(cfg.parcellation)==sel;
      fun(sel,:,:) = repmat(tmp(k,:,:), [sum(sel) 1]);
    end
  end
  source.(cfg.xparam) = source2.(cfg.xparam);
  if ~isempty(cfg.yparam)
    source.(cfg.yparam) = source2.(cfg.yparam);
  end
end

hasatlas = ~isempty(cfg.atlas);
if hasatlas
  if ischar(cfg.atlas)
    % initialize the atlas
    [p, f, x] = fileparts(cfg.atlas);
    fprintf(['reading ', f, ' atlas coordinates and labels\n']);
    atlas = ft_read_atlas(cfg.atlas);
  else
    atlas = cfg.atlas;
  end
end

hasparcels       = exist('parcels', 'var');
hasparcelslabel  = exist('parcelslabel', 'var');
hassurface       = isfield(cfg,'surffile') && ~isempty(cfg.surffile);

if size(source.pos)~=size(fun,1)
  error('inconsistent number of vertices in the cortical mesh');
end

if ~isfield(source, 'tri')
  if ~hassurface
    error('source.tri missing, this function requires a triangulated cortical sheet as source model');
  else% use surffile
      % read the triangulated cortical surface from file
      surf = ft_read_headshape(cfg.surffile);

      if isfield(surf, 'transform'),
        % compute the surface vertices in head coordinates
        surf.pos = ft_warp_apply(surf.transform, surf.pos);
      elseif ~isempty(cfg.surftrans),
        % compute the surface vertices in head coordinates
        surf.pos = ft_warp_apply(cfg.surftrans, surf.pos);
      end
      
      % convert units
      siz  = norm(idrange(source.pos));
      unit = ft_estimate_units(siz);
      surf = ft_convert_units(surf,unit);

      % downsample the cortical surface
      if cfg.surfdownsample > 1
        if ~isempty(cfg.surfinflated)
          error('downsampling the surface is not possible in combination with an inflated surface');
        end
        fprintf('downsampling surface from %d vertices\n', size(surf.pos,1));
        [temp.tri, temp.pos] = reducepatch(surf.tri, surf.pos, 1/cfg.surfdownsample);
        % find indices of retained patch faces
        [dummy, idx] = ismember(temp.pos, surf.pos, 'rows');
        idx(idx==0)  = [];
        surf.tri = temp.tri;
        surf.pos = temp.pos;
        clear temp
        % downsample other fields
        if isfield(surf, 'curv'),       surf.curv       = surf.curv(idx);       end
        if isfield(surf, 'sulc'),       surf.sulc       = surf.sulc(idx);       end
        if isfield(surf, 'hemisphere'), surf.hemisphere = surf.hemisphere(idx); end
      end

      % these are required
      if ~isfield(source, 'inside')
        source.inside = true(dim);
      end

      fprintf('%d voxels in functional data\n', prod(source.dim));
      fprintf('%d vertices in cortical surface\n', size(surf.pos,1));

      tmpcfg = [];
      tmpcfg.parameter = {cfg.funparameter};
      if ~isempty(cfg.maskparameter)
        tmpcfg.parameter = [tmpcfg.parameter {cfg.maskparameter}];
        maskparameter    = cfg.maskparameter;
      else
        tmpcfg.parameter = [tmpcfg.parameter {'mask'}];
        source.mask  = msk;
        maskparameter    = 'mask'; % temporary variable
      end
      tmpcfg.interpmethod = cfg.projmethod;
      tmpcfg.distmat      = cfg.distmat;
      tmpcfg.sphereradius = cfg.sphereradius;
      tmpcfg.projvec      = cfg.projvec;
      tmpcfg.projcomb     = cfg.projcomb;
      tmpcfg.projweight   = cfg.projweight;
      tmpcfg.projthresh   = cfg.projthresh;
      tmpdata             = ft_sourceinterpolate(tmpcfg, source, surf);

      if isfield(source, cfg.funparameter)
          fun = getsubfield(tmpdata, cfg.funparameter);
      end
      if issubfield(source, cfg.maskparameter)
          source.(cfg.maskparameter)    = tmpdata.(cfg.maskparameter);
      end
      
      % apply inflated
      if ~isempty(cfg.surfinflated)
       if ~isstruct(cfg.surfinflated)
         % read the inflated triangulated cortical surface from file
         surf = ft_read_headshape(cfg.surfinflated);
       else
         surf = cfg.surfinflated;
         if isfield(surf, 'transform'),
           % compute the surface vertices in head coordinates
           surf.pos = ft_warp_apply(surf.transform, surf.pos);
         elseif ~isempty(cfg.surftrans),
           % compute the surface vertices in head coordinates
           surf.pos = ft_warp_apply(cfg.surftrans, surf.pos);
         end
       end
      end
      
      % set surface data into source data
      source.pos    = surf.pos;
      source.tri    = surf.tri;
  end
end

if ~isempty(cfg.maskparameter) && ischar(cfg.maskparameter)
  mask = double(getsubfield(source, cfg.maskparameter));
else
  mask = 0.5*ones(size(fun));
end

xparam = source.(cfg.xparam);
if length(xparam)~=size(fun,2)
  error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.xparam);
end

if ~isempty(cfg.yparam)
  yparam = source.(cfg.yparam);
  if length(yparam)~=size(fun,3)
    error('inconsistent size of "%s" compared to "%s"', cfg.funparameter, cfg.yparam);
  end
else
  yparam = [];
end

if isempty(xlim)
  xlim(1) = min(xparam);
  xlim(2) = max(xparam);
end

xbeg = nearest(xparam, xlim(1));
xend = nearest(xparam, xlim(2));
% update the configuration
cfg.xlim = xparam([xbeg xend]);

if ~isempty(yparam)
  if isempty(ylim)
    ylim(1) = min(yparam);
    ylim(2) = max(yparam);
  end
  ybeg = nearest(yparam, ylim(1));
  yend = nearest(yparam, ylim(2));
  % update the configuration
  cfg.ylim = xparam([xbeg xend]);
  hasyparam = true;
else
  % this allows us not to worry about the yparam any more
  yparam = nan;
  ybeg = 1;
  yend = 1;
  cfg.ylim = [];
  hasyparam = false;
end

% make a subselection of the data
xparam  = xparam(xbeg:xend);
yparam  = yparam(ybeg:yend);
fun     = fun(:,xbeg:xend,ybeg:yend);
if hassource2 && isfield(source2, 'pos'),
  fun2 = fun2(:,xbeg:xend,ybeg:yend);
end
mask    = mask(:,xbeg:xend,ybeg:yend);
clear xbeg xend ybeg yend

if isempty(zlim)
  zlim(1) = min(fun(:));
  zlim(2) = max(fun(:));
  if zlim(1)>=zlim(2)
    zlim(1) = 0;
    zlim(2) = 1;
  end
  % set 0 for too small value
  if abs(zlim(1)) <= abs(zlim(2)) * 10^-5
      zlim(1) = 0;
  elseif abs(zlim(2)) <= abs(zlim(1)) * 10^-5
      zlim(2) = 0;
  end
  % update the configuration
  cfg.zlim = zlim;
end

if isempty(olim)
  olim(1) = min(mask(:));
  olim(2) = max(mask(:));
  if olim(1)>=olim(2)
    olim(1) = 0;
    olim(2) = 1;
  end
  % set 0 for too small value
  if abs(olim(1)) <= abs(olim(2)) * 10^-5
      olim(1) = 0;
  elseif abs(olim(2)) <= abs(olim(1)) * 10^-5
      olim(2) = 0;
  end
  % update the configuration
  cfg.alim = olim;
end

% collect the data and the options to be used in the figure
opt.cfg     = cfg;
opt.xparam  = xparam;
opt.yparam  = yparam;
opt.xval    = 0;
opt.yval    = 0;
opt.dat     = fun;
opt.mask    = abs(mask);
opt.pos     = source.pos;
opt.tri     = source.tri;
if isfield(source, 'inside')
  opt.vindx   = source.inside(:);
else
  opt.vindx   = 1:size(opt.pos,1);
end
opt.speed   = 1;
opt.record  = 0;
opt.threshold = 0;
opt.frame   = 0;
opt.cleanup = false;
if hasparcels,      opt.parcellation      = parcels; end
if hasparcelslabel, opt.parcellationlabel = parcelslabel; end
if hasatlas,
    opt.atlas       = atlas;
    opt.queryrange  = cfg.queryrange;
    opt.coordsys    = source.coordsys;
end

% add functional data of optional third input to the opt structure
% FIXME here we should first check whether the meshes correspond!
if hassource2 && isfield(source2, 'pos')
  opt.dat2 = fun2;
  opt.dat1 = opt.dat;
end

%% start building the figure
h = figure;
set(h, 'color', [1 1 1]);
set(h, 'visible', 'on');
set(h, 'renderer', cfg.renderer);
set(h, 'toolbar', 'figure');
set(h, 'CloseRequestFcn', @cb_quitbutton);
set(h, 'position', [100 200 700 500]);
set(h, 'WindowButtonDownFcn', @cb_getposition);
%-- release 3D rotation
set(h, 'WindowButtonUpFcn', @cb_3drotfin);
set(h, 'WindowButtonMotionFcn', @cb_3drotupdate);
if ~isempty(cfg.title)
  title(cfg.title);
end

% get timer object
t = timer;
set(t, 'timerfcn', {@cb_timer, h}, 'period', 0.1, 'executionmode', 'fixedSpacing');

% make the user interface elements
cambutton    = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'light', 'userdata', 'C');
playbutton   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'play',   'userdata', 'p');
recordbutton = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'record', 'userdata', 'r');
quitbutton   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'quit',   'userdata', 'q');
if isfield(opt, 'dat2'),
  displaybutton = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'display: var1',   'userdata', 'f');
end

thrmin   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'downarrow');
thr      = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'threshold', 'userdata', 't');
thrplus  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'uparrow');
spdmin   = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'shift+downarrow');
spd      = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'speed','userdata', 's');
spdplus  = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'shift+uparrow');
olim       = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'masklim', 'userdata', 'a');
olimminmin = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'ctrl+leftarrow');
olimmaxmin = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'ctrl+shift+leftarrow');
olimminplus = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'ctrl+rightarrow');
olimmaxplus = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'ctrl+shift+rightarrow');
clim       = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'colorlim', 'userdata', 'z');
climminmin = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'leftarrow');
climmaxmin = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+leftarrow');
climminplus = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'rightarrow');
climmaxplus = uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+rightarrow');
sliderx  = uicontrol('parent', h, 'units', 'normalized', 'style', 'slider',     'string', sprintf('%s = ', cfg.xparam));
stringx  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
slidery  = uicontrol('parent', h, 'units', 'normalized', 'style', 'slider',     'string', sprintf('%s = ', cfg.yparam));
stringy  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
stringz  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
stringp  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');
stringa  = uicontrol('parent', h, 'units', 'normalized', 'style', 'text');

if isfield(opt,'dat2')
  set(displaybutton, 'position', [0.005 0.34 0.18 0.05], 'callback', @cb_keyboard);
end

set(cambutton,    'position', [0.095 0.34 0.09 0.05], 'callback', @cb_keyboard);
set(quitbutton,   'position', [0.005 0.34 0.09 0.05], 'callback', @cb_keyboard);
set(playbutton,   'position', [0.005 0.28 0.09 0.05], 'callback', @cb_keyboard);
set(recordbutton, 'position', [0.095 0.28 0.09 0.05], 'callback', @cb_keyboard);
set(thrmin,       'position', [0.005 0.22 0.03 0.05], 'callback', @cb_keyboard);
set(thr,          'position', [0.035 0.22 0.12 0.05], 'callback', @cb_keyboard);
set(thrplus,      'position', [0.155 0.22 0.03 0.05], 'callback', @cb_keyboard);
set(olimminmin,   'position', [0.005 0.16  0.03 0.025], 'callback', @cb_keyboard);
set(olimmaxmin,   'position', [0.005 0.185 0.03 0.025], 'callback', @cb_keyboard);
set(olim,         'position', [0.035 0.16 0.12 0.05], 'callback', @cb_keyboard);
set(olimminplus,  'position', [0.155 0.16  0.03 0.025], 'callback', @cb_keyboard);
set(olimmaxplus,  'position', [0.155 0.185 0.03 0.025], 'callback', @cb_keyboard);
set(climminmin,   'position', [0.005 0.10  0.03 0.025], 'callback', @cb_keyboard);
set(climmaxmin,   'position', [0.005 0.125 0.03 0.025], 'callback', @cb_keyboard);
set(clim,         'position', [0.035 0.10 0.12 0.05], 'callback', @cb_keyboard);
set(climminplus,  'position', [0.155 0.10  0.03 0.025], 'callback', @cb_keyboard);
set(climmaxplus,  'position', [0.155 0.125 0.03 0.025], 'callback', @cb_keyboard);
set(spdmin,       'position', [0.005 0.04 0.03 0.05], 'callback', @cb_keyboard);
set(spd,          'position', [0.035 0.04 0.12 0.05], 'callback', @cb_keyboard);
set(spdplus,      'position', [0.155 0.04 0.03 0.05], 'callback', @cb_keyboard);
set(sliderx,      'position', [0.02 0.45 0.3 0.03], 'callback',  @cb_slider);%[0.200 0.04  0.78 0.03], 'callback', @cb_slider);
set(slidery,      'position', [0.350 0.55  0.03 0.35], 'callback', @cb_slider);
set(stringx,      'position', [0.800 0.93 0.18 0.03]);
set(stringy,      'position', [0.800 0.90 0.18 0.03]);
set(stringz,      'position', [0.650 0.96 0.33 0.03]);
set(stringp,      'position', [0.650 0.87 0.33 0.03]);
set(stringa,      'position', [0.700 0.90-0.03*(hasyparam+hasparcels) 0.28 0.03]);      % set based on the exsitence of parcels

set(stringx, 'string', sprintf('%s = ', cfg.xparam));
set(stringy, 'string', sprintf('%s = ', cfg.yparam));
set(stringz, 'string', sprintf('position = '));
set(stringp, 'string', sprintf('parcel = '));
set(stringa, 'string', sprintf('atlas = '));
set(stringx, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringy, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringz, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringp, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);
set(stringa, 'horizontalalignment', 'right', 'backgroundcolor', [1 1 1]);

% create axes object to contain the mesh
hx = axes;
set(hx, 'position', [0.4 0.08 0.6 0.8]);
set(hx, 'tag', 'mesh');
if isfield(source, 'sulc')
  vdat = source.sulc;
  vdat = vdat-min(vdat);
  vdat = vdat./max(vdat);
  vdat = 0.1+0.3.*repmat(round(1-vdat),[1 3]);
  hs1 = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', vdat);
else
  hs1 = ft_plot_mesh(source, 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]);
end
lighting gouraud
siz = [size(opt.dat) 1];
hs = ft_plot_mesh(source, 'edgecolor', 'none', 'vertexcolor', 0*opt.dat(:,ceil(siz(2)/2),ceil(siz(3)/2)), 'facealpha', 0*opt.mask(:,1,1));
lighting gouraud
cam1 = camlight('left');
cam2 = camlight('right');
caxis(cfg.zlim);
alim(cfg.alim);
colormap(cfg.colormap);

% create axis object to contain a time course
hy = axes;
set(hy, 'position', [0.02 0.55 0.3 0.35]);
set(hy, 'yaxislocation', 'right');

if ~hasyparam
  tline = plot(opt.xparam, mean(opt.dat(opt.vindx,:))); hold on;
  abc = axis;
  axis([opt.xparam(1) opt.xparam(end) abc(3:4)]);
  vline = plot(opt.xparam(1)*[1 1], abc(3:4), 'r');

  if hassource2 && isfield(source2, 'pos')
    tline2 = plot(opt.xparam, mean(opt.dat2(opt.vindx,:)), 'r'); hold on;
  end

else
  tline = imagesc(opt.xparam, opt.yparam, shiftdim(mean(opt.dat(opt.vindx,:,:)),1)'); axis xy; hold on;
  abc   = [opt.xparam([1 end]) opt.yparam([1 end])];
  vline = plot(opt.xparam(ceil(siz(2)/2)).*[1 1], abc(3:4));
  hline = plot(abc(1:2), opt.yparam(ceil(siz(3)/2)).*[1 1]);
  %error('not yet implemented');
end
set(hy, 'tag', 'timecourse');

% remember the various handles
opt.h   = h;  % handle to the figure
opt.hs  = hs; % handle to the mesh
opt.hx  = hx; % handle to the axes containing the mesh
opt.hy  = hy; % handle to the axes containing the timecourse
opt.cam = [cam1 cam2]; % handles to the light objects
opt.vline = vline; % handle to the line in the ERF plot
opt.tline = tline; % handle to the ERF
if exist('hline', 'var')
  opt.hline = hline;
end
if hassource2 && isfield(source2, 'pos'),
  opt.tline2 = tline2;
end
opt.playbutton   = playbutton; % handle to the playbutton
opt.recordbutton = recordbutton; % handle to the recordbutton
opt.quitbutton   = quitbutton; % handle to the quitbutton
try, opt.displaybutton = displaybutton; end

%opt.p   = p;
opt.t   = t;
%opt.hx  = hx;
%opt.hy  = hy;
opt.sliderx  = sliderx;
opt.slidery  = slidery;
opt.stringx  = stringx;
opt.stringy  = stringy;
opt.stringz  = stringz;
opt.stringp  = stringp;
opt.stringa  = stringa;

if ~hasyparam
  set(opt.slidery, 'visible', 'off');
  set(opt.stringy, 'visible', 'off');
end

if ~hasparcels
  set(opt.stringp, 'visible', 'off');
end
if ~hasatlas
  set(opt.stringa, 'visible', 'off');
end

opt.nargout = nargout;
setappdata(h, 'opt', opt);
cb_slider(h);     set(h,'Pointer','arrow');

if nargout
    while opt.cleanup==0
      uiwait(h);
      opt = getappdata(h, 'opt');
    end
    stop(opt.t);

    if isfield(opt,'movie')
      M = opt.movie;
    else
      M = [];
    end

    delete(h);

end

ft_private('reset','silent');

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous source
ft_postamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_slider(h, eventdata)

persistent previous_valx previous_valy previous_vindx previous_threshold

if isempty(previous_valx)
  previous_valx = 0;
end
if isempty(previous_valy)
  previous_valy = 0;
end

h    = getparent(h);
opt  = getappdata(h, 'opt');
valx = get(opt.sliderx, 'value');
valx = round(valx*(size(opt.dat,2)-1))+1;
valx = min(valx, size(opt.dat,2));
valx = max(valx, 1);

valy = get(opt.slidery, 'value');
valy = round(valy*(size(opt.dat,3)-1))+1;
valy = min(valy, size(opt.dat,3));
valy = max(valy, 1);

mask = opt.mask(:,valx,valy);
mask(opt.dat(:,valx,valy)<opt.threshold) = 0;

% update stuff
if previous_valx~=valx || previous_valy~=valy || previous_threshold ~= opt.threshold
  % update strings
  set(opt.stringx, 'string', sprintf('%s = %3.3f\n', opt.cfg.xparam, opt.xparam(valx)));
  set(opt.stringy, 'string', sprintf('%s = %3.3f\n', opt.cfg.yparam, opt.yparam(valy)));

  % update data in mesh
  set(opt.hs, 'FaceVertexCData',     opt.dat(:,valx,valy));
  set(opt.hs, 'FaceVertexAlphaData', mask);

  set(opt.vline, 'xdata', [1 1]*opt.xparam(valx));
  if isfield(opt, 'hline')
    set(opt.hline, 'ydata', [1 1]*opt.yparam(valy));
  end
end

% update ERF-plot
if ~isfield(opt, 'hline')
  set(opt.hy,    'ylim',   opt.cfg.zlim);
  set(opt.vline, 'ydata',  opt.cfg.zlim);
else
  set(opt.hy,    'clim',   opt.cfg.zlim);
end
if ~(numel(previous_vindx)==numel(opt.vindx) && all(previous_vindx==opt.vindx))
  if ~isfield(opt, 'hline')
    tmp = mean(opt.dat(opt.vindx,:,valy),1);
    set(opt.tline, 'ydata', tmp);
  else
    tmp = shiftdim(mean(opt.dat(opt.vindx,:,:),1))';
    set(opt.tline, 'cdata', tmp);
  end
  %set(opt.hy,    'ylim',  [min(tmp(:)) max(tmp(:))]);
  %set(opt.vline, 'ydata', [min(tmp(:)) max(tmp(:))]);

  if isfield(opt, 'dat2')
    tmp = mean(opt.dat1(opt.vindx,:,valy),1);
    set(opt.tline, 'ydata', tmp);
    tmp = mean(opt.dat2(opt.vindx,:,valy),1);
    set(opt.tline2, 'ydata', tmp);
  end

  set(opt.hy,    'yaxislocation', 'right');
  set(opt.stringz, 'string', sprintf('position = [%2.1f, %2.1f, %2.1f]', opt.pos(opt.vindx,:)));
  if isfield(opt, 'parcellation'),
    set(opt.stringp, 'string', sprintf('parcel = %s', opt.parcellationlabel{opt.parcellation(opt.vindx)}));
  end
  if isfield(opt, 'atlas'),
    lab = load_atlas(opt);
    set(opt.stringa, 'string', sprintf('atlas = %s', lab));
  end
end

if opt.record
  tmp = get(opt.h, 'position');
  opt.frame = opt.frame + 1;
  opt.movie(opt.frame) = getframe(opt.h,[1 1 tmp(3:4)-1]);
end
setappdata(h, 'opt', opt);

previous_valx = valx;
previous_valy = valy;
previous_vindx = opt.vindx;
previous_threshold = opt.threshold;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_playbutton(h, eventdata)

opt = getappdata(h, 'opt');
if strcmp(get(opt.playbutton, 'string'), 'pause')
  stop(opt.t);
  set(opt.playbutton, 'string', 'play');
else
  start(opt.t);
  set(opt.playbutton, 'string', 'pause');
end
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quitbutton(h, eventdata)

opt = getappdata(h, 'opt');
opt.cleanup = 1;
setappdata(h, 'opt', opt);
if strcmp(get(opt.playbutton, 'string'), 'pause')
    cb_playbutton(h);   % stop playing
end
if strcmp(get(opt.recordbutton, 'string'), 'stop')
    cb_recordbutton(h);   % stop recording
end
if ~opt.nargout
    h   = getparent(h);
    delete(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_recordbutton(h, eventdata)

opt = getappdata(h, 'opt');
if strcmp(get(opt.recordbutton, 'string'), 'stop')
  opt.record = 0;
  set(opt.recordbutton, 'string', 'record');
else
  opt.record = 1;
  set(opt.recordbutton, 'string', 'stop');
end
setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_timer(obj, info, h)

opt   = getappdata(h, 'opt');
delta = opt.speed/size(opt.dat,2);
val = get(opt.sliderx, 'value');
val = val + delta;
if val>1
  val = val-1;
end
set(opt.sliderx, 'value', val);
setappdata(h, 'opt', opt);
cb_slider(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_alim(h, eventdata)
if ~ishandle(h)
  return
end
opt = guidata(h);
switch get(h, 'String')
  case '+'
    alim(alim*sqrt(2));
  case '-'
    alim(alim/sqrt(2));
end % switch
guidata(h, opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)
ft_private('silent');

h   = getparent(h);
opt = getappdata(h, 'opt');
lastmousebttn = get(h,'selectiontype');
if strcmp(get(get(h, 'currentaxes'), 'tag'), 'timecourse')
  % get the current point
  pos = get(opt.hy, 'currentpoint');
  set(opt.sliderx, 'value', nearest(opt.xparam, pos(1,1))./numel(opt.xparam));
  if isfield(opt, 'hline')
    set(opt.slidery, 'value', nearest(opt.yparam, pos(1,2))./numel(opt.yparam));
  end
elseif strcmp(get(get(h, 'currentaxes'), 'tag'), 'mesh') 
  switch lastmousebttn
      case 'alt'
      % 3D rotation for right-click & drag
      caid = get(h, 'currentaxes');
      
      opt.holding                            = true;
      opt.previous_pos                       = get(groot,'PointerLocation');
      [opt.viewpoint(1), opt.viewpoint(2)]   = view(caid);
      
      cur_unit = get(caid,'Units');
      set(caid,'Units','pixels');
      gco_res  = get(caid,'Position');
      set(caid,'Units',cur_unit);
      
      opt.mvspeed = 180 ./ min(gco_res(3:4));
      
      set(h,'Pointer','fleur');
      
      otherwise
      if ~isfield(opt,'holding') || ~opt.holding
      % get the current point, which is defined as the intersection through the
      % axis-box (in 3D)
      pos       = get(opt.hx, 'currentpoint');

      % get the intersection with the mesh
      [ipos, d] = intersect_line(opt.pos, opt.tri, pos(1,:), pos(2,:));
      if ~isempty(ipos)
          [md, ix]  = min(abs(d));

          dpos      = opt.pos - ipos(ix*ones(size(opt.pos,1),1),:);
          opt.vindx = nearest(sum(dpos.^2,2),0);
      end
      end
  end
end
setappdata(h, 'opt', opt);
ft_private('reset','silent');
cb_slider(h);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(h,eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h = getparent(h);
opt = getappdata(h, 'opt');

switch key
  case 'leftarrow' % change colorlim
    opt.cfg.zlim(1) = opt.cfg.zlim(1)-0.1*abs(opt.cfg.zlim(1));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.zlim);
    set(opt.hx, 'Clim', opt.cfg.zlim);

  case 'shift+leftarrow' % change colorlim
    opt.cfg.zlim(1) = opt.cfg.zlim(1)+0.1*abs(opt.cfg.zlim(1));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.zlim);
    set(opt.hx, 'Clim', opt.cfg.zlim);

  case 'rightarrow'
    opt.cfg.zlim(2) = opt.cfg.zlim(2)-0.1*abs(opt.cfg.zlim(2));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.zlim);
    set(opt.hx, 'Clim', opt.cfg.zlim);

  case 'shift+rightarrow'
    opt.cfg.zlim(2) = opt.cfg.zlim(2)+0.1*abs(opt.cfg.zlim(2));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.zlim);
    set(opt.hx, 'Clim', opt.cfg.zlim);

  case 'ctrl+leftarrow' % change opacitylim
    opt.cfg.alim(1) = opt.cfg.alim(1)-0.1*abs(opt.cfg.alim(1));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.alim);
    set(opt.hx, 'Alim', opt.cfg.alim);

  case 'ctrl+shift+leftarrow' % change opacitylim
    opt.cfg.alim(1) = opt.cfg.alim(1)+0.1*abs(opt.cfg.alim(1));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.alim);
    set(opt.hx, 'Alim', opt.cfg.alim);

  case 'ctrl+rightarrow'
    opt.cfg.alim(2) = opt.cfg.alim(2)-0.1*abs(opt.cfg.alim(2));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.alim);
    set(opt.hx, 'Alim', opt.cfg.alim);

  case 'ctrl+shift+rightarrow'
    opt.cfg.alim(2) = opt.cfg.alim(2)+0.1*abs(opt.cfg.alim(2));
    setappdata(h, 'opt', opt);
    caxis(opt.cfg.alim);
    set(opt.hx, 'Alim', opt.cfg.alim);

  case 'uparrow' % enhance threshold
    opt.threshold = opt.threshold+0.01.*max(opt.dat(:));
    setappdata(h, 'opt', opt);
  case 'downarrow' % lower threshold
    opt.threshold = opt.threshold-0.01.*max(opt.dat(:));
    setappdata(h, 'opt', opt);
  case 'shift+uparrow' % change speed
    opt.speed = opt.speed*sqrt(2);
    setappdata(h, 'opt', opt);
  case 'shift+downarrow'
    opt.speed = opt.speed/sqrt(2);
    opt.speed = max(opt.speed, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
  case 'ctrl+uparrow' % change channel
  case 'C' % update camera position
    camlight(opt.cam(1), 'left');
    camlight(opt.cam(2), 'right');
  case 'p'
    cb_playbutton(h);
  case 'q'
    cb_quitbutton(h);
  case 'r'
    cb_recordbutton(h);
  case 's'
    % select the speed
    response = inputdlg('speed', 'specify', 1, {num2str(opt.speed)});
    if ~isempty(response)
      opt.speed = str2double(response);
      setappdata(h, 'opt', opt);
    end
  case 't'
    % select the threshold
    response = inputdlg('threshold', 'specify', 1, {num2str(opt.threshold)});
    if ~isempty(response)
      opt.threshold = str2double(response);
      setappdata(h, 'opt', opt);
    end
  case 'z'
    % select the colorlim
    response = inputdlg('colorlim', 'specify', 1, {[num2str(opt.cfg.zlim(1)),' ',num2str(opt.cfg.zlim(2))]});
    if ~isempty(response)
      [tok1, tok2] = strtok(response, ' ');
      opt.cfg.zlim(1) = str2double(deblank(tok1));
      opt.cfg.zlim(2) = str2double(deblank(tok2));
      set(opt.hx, 'Clim', opt.cfg.zlim);
      setappdata(h, 'opt', opt);
    end
  case 'a'
    % select the opacitylim
    response = inputdlg('opacitylim for mask', 'specify', 1, {[num2str(opt.cfg.alim(1)),' ',num2str(opt.cfg.alim(2))]});
    if ~isempty(response)
      [tok1, tok2] = strtok(response, ' ');
      opt.cfg.alim(1) = str2double(deblank(tok1));
      opt.cfg.alim(2) = str2double(deblank(tok2));
      set(opt.hx, 'Alim', opt.cfg.alim);
      setappdata(h, 'opt', opt);
    end
  case 'f'
    if isfield(opt, 'dat2')
      if isequaln(opt.dat,opt.dat2),
        opt.dat = opt.dat1;
        set(opt.displaybutton, 'string', 'display: var1');
      end
      if isequaln(opt.dat,opt.dat1),
        opt.dat = opt.dat2;
        set(opt.displaybutton, 'string', 'display: var2');
      end
    end
    setappdata(h, 'opt', opt);
    cb_slider(h);
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
%     cb_help(h);
end

if ishandle(h)
    cb_slider(h);
    uiresume(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(handle, eventdata)
ft_private('silent');
if (isempty(eventdata) && ft_platform_supports('matlabversion',-Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  %-- determine the key that corresponds to the uicontrol element that was activated
  key = get(handle, 'userdata');
else
  %-- determine the key that was pressed on the keyboard
  key = eventdata.Key;
  
    % handle possible numpad events (different for Windows and UNIX systems)
    % NOTE: shift+numpad number does not work on UNIX, since the shift
    % modifier is always sent for numpad events
    if isunix()
      shiftInd = match_str(eventdata.Modifier, 'shift');
      if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
        % now we now it was a numpad keystroke (numeric character sent AND
        % shift modifier present)
        key = eventdata.Character;
        eventdata.Modifier(shiftInd) = []; % strip the shift modifier
      end
    elseif ispc()
      if strfind(eventdata.Key, 'numpad')
        key = eventdata.Character;
      end
    end

    if ~isempty(eventdata.Modifier)
      key = [eventdata.Modifier{1} '+' key];
    end

end
ft_private('reset','silent');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_3drotfin(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');
lastmousebttn = get(h,'selectiontype');

if isfield(opt,'holding') && opt.holding
    switch lastmousebttn
        case 'alt'
            opt.holding = false;
    end
    set(h,'Pointer','arrow');
end
setappdata(h, 'opt', opt);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_3drotupdate(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

if isfield(opt,'holding') && opt.holding
    pos       = get(groot,'PointerLocation');

    moving    = opt.previous_pos - pos;
    moving    = mean(moving(:,1:2),1);
        
    newVP     = [opt.viewpoint + moving .* opt.mvspeed];
    newVP(2)  = max(min(newVP(2),90),-90);
    
    view(get(h, 'currentaxes'),newVP);
    
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lab = load_atlas(opt)

ft_private('silent');
if length(opt.vindx)==1
    % determine the anatomical label of the current position
    lab = atlas_lookup(opt.atlas, opt.pos(opt.vindx,:), 'inputcoord', opt.coordsys, 'queryrange', opt.queryrange);
    if isempty(lab)
        lab = 'NA';
    else
        tmp = sprintf('%s', strrep(lab{1}, '_', ' '));
        for i=2:length(lab)
            tmp = [tmp sprintf(', %s', strrep(lab{i}, '_', ' '))];
        end
        lab = tmp;
    end
else
    lab = 'NA';
end
ft_private('reset','silent');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDRANGE  range for more robust range estimation (copied from ft_convert_units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = idrange(x)
keeprow=true(size(x,1),1);
for l=1:size(x,2)
  keeprow = keeprow & isfinite(x(:,l));
end
sx = sort(x(keeprow,:), 1);
ii = round(interp1([0, 1], [1, size(x(keeprow,:), 1)], [.1, .9]));  % indices for 10 & 90 percentile
r = diff(sx(ii, :));
