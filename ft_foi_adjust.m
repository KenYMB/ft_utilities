function [cfg, fbin, padding, data] = ft_foi_adjust(cfg, data)

% FT_FOI_ADJUST adjusts cfg.foi values considering padding for fourier transform.
% If cfg.foilim is specified instead of cfg.foi, it calculates correct foi. 
% If you want to apply padding automatically on data, 'data' should be returned.
% 
% Use as
%   [cfg] = ft_foi_adjust(cfg)
%   [cfg] = ft_foi_adjust(cfg, data)
%   [cfg, fbin, padding] = ft_foi_adjust(cfg, data)
%   [cfg, fbin, padding, data] = ft_foi_adjust(cfg, data)
% 
% where
%   fbin        is an applyed frequency bin.
%   padding     is padded data length.
% 
% The following configuration options are needed:
%   cfg.foi        = vector 1 x numfoi, frequencies of interest
%       OR
%   cfg.foilim     = [begin end], frequency band of interest
% 
% The following configuration options are needed if 'data' does not exist:
%   cfg.fsample    = number, a sampling frequency
%   cfg.trllength  = number, a length of time points
%                    you can specify cfg.padding instead of cfg.trllength
% These values are prior to 'data'
% 
% The following configuration options are supported:
%   cfg.pad        = number, 'nextpow2' (default) 'maxperlen', or 'nextinteger'
%                       'nextinteger': set padded data length as integer in second
%                                      frequency steps come to integer
%   cfg.t_ftimwin  = vector 1 x numfoi, length of time window (in seconds)
%   
%   cfg.correctt_ftimwin = 'yes'(default) or 'no'
%                          make new 'cfg.t_ftwin' based on input cfg.t_ftimwin and updated cfg.foi
%   cfg.wave_ftimwin  = number, wave number for time window (available when cfg.correctt_ftimwin is 'yes')
%                       this option have priority to cfg.t_ftwin
% 
% If the following configuration option is input, cfg.pad is calculated and replaced: 
%   cfg.padding    = number, 'maxperlen', or 'maxftimwin'
%                    length (in seconds) to which the trials are padded for frequency analysis (default = 0)
%                    padding length is calculated based on cfg.padding instead of data length
%                    Example: padding length is calculated based on data length of 2s
%                           cfg.pad         = 'nextpow2';
%                           cfg.padding     = 2;
%                           cfg             = ft_foi_adjust(cfg);
%                    Example: calculate padding and apply on data
%                           cfg.pad         = 'nextpow2';
%                           cfg.t_ftimwin   = 5./cfg.foi;
%                           cfg.padding     = cfg.t_ftimwin(1) + diff(cfg.toi([1,end]));
%                                             % same as 'maxftimwin'
%                           [cfg,~,~,data]  = ft_foi_adjust(cfg,data);
% 
% The following configuration options are supported when 'data' is specified as a return value
%   cfg.padtype    = 'zero', 'mean', 'localmean', 'edge', 'mirror', 'nan' or 'remove' (default = 'zero')
%                    if {pre_padtype, post_padtype} is specified, pre_pdtype is returned as cfg.padtype
%   cfg.paddir     = direction of padding, 'left'/'right'/'both' (default = 'both')
% 
% See also FT_APPLY_PADDING, FT_FREQANALYSIS

% using: fieldtrip, ft_apply_padding

% 20161220 Yuasa: create for ft_freqanalysis
% 20161221 Yuasa: minor fix
% 20170217 Yuasa: remove repeted foi & enable padding length
% 20170220 Yuasa: add ft_apply_padding
% 20170221 Yuasa: change usage of cfg.padding
% 20170406 Yuasa: enable correctt_ft_win for foilim input
% 20170410 Yuasa: check if 'mtmconvol' for padding-apply
% 20180424 Yuasa: add 'cfg.wave_ftimwin', add 'maxftimwin' for 'cfg.padding'
%                 add 'nextinteger' for 'cfg.pad'

narginchk(1,2);

ft_defaults;

% set all the defaults
cfg.pad       = ft_getopt(cfg, 'pad',       []);
if isempty(cfg.pad)
  warning('Default cfg.pad = ''nextpow2'' is set.  It is different from the defalt value of ft_freqanalysis.');
  cfg.pad = 'nextpow2';
end
cfg.foi       = ft_getopt(cfg, 'foi',       []);
cfg.foilim    = ft_getopt(cfg, 'foilim',    []);
cfg.correctt_ftimwin = ft_getopt(cfg, 'correctt_ftimwin', 'yes');
cfg.padding = ft_getopt(cfg, 'padding', 'maxperlen');

%-- if data exist
if nargin == 2
  %-- check if the input data is valid for this function
  data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'mvar'});
  %-- determine sampling frequency
  fsample = data.fsample;
  %-- determine trial characteristics
  trllength = zeros(1, numel(data.trial));
  for itrial = 1:numel(data.trial)
     trllength(itrial) = size(data.trial{itrial}, 2);
  end
end

%-- set parameter
if isfield(cfg,'fsample') && ~isempty(cfg.fsample)
  fsample = cfg.fsample;
end
assert(logical(exist('fsample','var')),'you must specify a sampling frequency');
if isfield(cfg,'trllength') && ~isempty(cfg.trllength)
  if exist('trllength','var'), warning('data length is replaced by cfg.trllength'); end
  trllength = cfg.trllength;
end
if isnumeric(cfg.padding) && ~isempty(cfg.padding)
  padlen    = cfg.padding * fsample;
  if exist('trllength','var') && padlen < max(trllength)
    warning('the specified padding (cfg.padding) is too short');
  else    
    trllength = padlen;
  end
elseif strcmp(cfg.padding, 'maxftimwin')
    if ~isempty(cfg.foilim),    fmin = min(cfg.foilim(:));
    else                        fmin = min(cfg.foi(:));
    end
    padlen = 0;
    if isfield(cfg,'wave_ftimwin') && ~isempty(cfg.wave_ftimwin)
      padlen = max(cfg.wave_ftimwin) ./ fmin .* fsample;
    elseif isfield(cfg,'t_ftimwin') && ~isempty(cfg.t_ftimwin)
      padlen = max(cfg.t_ftimwin(:)) .* fsample;
    end
    trllength = trllength + padlen;
elseif ~strcmp(cfg.padding, 'maxperlen')
    error('cfg.padding is invalid');
end
assert(logical(exist('trllength','var')),'you must specify a length of time points');

%-- padding
if strcmp(cfg.pad, 'maxperlen')
  padding = max(trllength);
  pad = padding/fsample;
elseif strcmp(cfg.pad, 'nextpow2')
  padding = 2^nextpow2(max(trllength));
  pad = padding/fsample;
elseif strcmp(cfg.pad, 'nextinteger')
  padding = max(trllength);
  pad = ceil(padding/fsample);
  padding = pad*fsample;
else
  warning('cfg.pad is numeric. cfg.padding is ignored.');
  pad = cfg.pad;
  padding = pad*fsample;
  if padding<max(trllength)
    error('the specified padding (cfg.pad) is too short');
  end
end

%-- new foi
if ~isempty(cfg.foi) && ~isempty(cfg.foilim)
  error('use either cfg.foi or cfg.foilim')
elseif ~isempty(cfg.foilim)
  %-- get the full foi in the current foilim and set it too be used as foilim
  fboilim = round(cfg.foilim .* pad) + 1;
  fboi    = fboilim(1):1:fboilim(2);
  cfg.foi = (fboi-1) ./ pad;
  fprintf(1,'cfg.foi is now adjusted.\n');
  
  if strcmp(cfg.correctt_ftimwin,'yes')
    if isfield(cfg,'wave_ftimwin') && ~isempty(cfg.wave_ftimwin)
      cfg.t_ftimwin = cfg.wave_ftimwin(1) ./ cfg.foi;
      fprintf(1,'cfg.t_ftimwin is now prepared.\n');
    elseif isfield(cfg,'t_ftimwin') && ~isempty(cfg.t_ftimwin)
      cfg.t_ftimwin = cfg.t_ftimwin(1) * ones(1,length(cfg.foi));
      fprintf(1,'cfg.t_ftimwin is now adjusted.\n');
    else
      warning('t_ftimwin correction is skipped.  cfg.t_ftimwin is not specified.');
    end
  end
else
  %-- correct foi if foilim was empty and try to correct t_ftimwin (by detecting whether there is a constant factor between foi and t_ftimwin: cyclenum)
  oldfoi = cfg.foi;
  [fboi, unia, unib]   = unique(round(cfg.foi .* pad) + 1); % remove repeated foi
  repidx               = setdiff(1:length(unib),unia);      % specify removed index
  cfg.foi    = (fboi-1) ./ pad; % boi - 1 because 0 Hz is included in fourier output
  fprintf(1,'cfg.foi is now adjusted.\n');
  
  if strcmp(cfg.correctt_ftimwin,'yes')
    if isfield(cfg,'wave_ftimwin') && ~isempty(cfg.wave_ftimwin)
      cfg.t_ftimwin = cfg.wave_ftimwin(1) ./ cfg.foi;
      fprintf(1,'cfg.t_ftimwin is now prepared.\n');
    elseif isfield(cfg,'t_ftimwin') && ~isempty(cfg.t_ftimwin)
      if all(cfg.t_ftimwin==cfg.t_ftimwin(1));
         cfg.t_ftimwin = cfg.t_ftimwin(1) * ones(1,length(cfg.foi));
         fprintf(1,'cfg.t_ftimwin is now adjusted.\n');
      elseif length(oldfoi) == length(cfg.t_ftimwin)
        cyclenum = oldfoi .* cfg.t_ftimwin;
        if all((cyclenum - cyclenum(1)) < (20*eps))
          cyclenum(repidx) = [];                                % remove repeated foi
          cfg.t_ftimwin = cyclenum ./ cfg.foi;
          fprintf(1,'cfg.t_ftimwin is now adjusted.\n');
        else
          warning('t_ftimwin correction is skipped.  failed to identify t_ftimwin.');
        end
      else
        warning('t_ftimwin correction is skipped.  failed to identify t_ftimwin.');
      end
    else
      warning('t_ftimwin correction is skipped.  cfg.t_ftimwin is not specified.');
    end
  end
end
cfg     = rmfield(cfg,'foilim');
cfg     = rmfield(cfg,'correctt_ftimwin');
fbin    = 1 / pad;


%-- new pad
if isfield(cfg,'padding') && ~isempty(cfg.padding)
    cfg.pad      = pad;
end

%-- apply padding
if nargout >= 4
  nomethod = ~isfield(cfg,'method') || isempty(cfg.method);
  if nomethod || strcmp(cfg.method,'mtmconvol') || strcmp(cfg.method,'wavelet') || strcmp(cfg.method,'tfr') || strcmp(cfg.method,'hilbert')
    if exist('ft_apply_padding')==2
      tmpcfg          = [];
      tmpcfg.padding  = pad;
      tmpcfg.padtype  = ft_getopt(cfg, 'padtype',      []);
      tmpcfg.paddir   = ft_getopt(cfg, 'paddir',       []);
      data            = ft_apply_padding(tmpcfg,data);
      fprintf(1+nomethod,'padding is applied.\n');
    else
      warning('''ft_apply_padding'' is not available. Skip applying padding.');
    end
  else
    warning('Padding in advance is not recommended for ''%s''.  Skip applying padding.',cfg.method);
  end
end

%-- cfg reset
if isfield(cfg,'padding')
    cfg          = rmfield(cfg,'padding');
end
if isfield(cfg,'paddir')
    cfg          = rmfield(cfg,'paddir');
end
if isfield(cfg,'padtype') && iscell(cfg.padtype)
    cfg.padtype  = cfg.padtype{1};
end
if isfield(cfg,'wave_ftimwin')
    cfg          = rmfield(cfg,'wave_ftimwin');
end
