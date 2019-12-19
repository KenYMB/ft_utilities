function [data] = ft_interpTFR(cfg, data)

% FT_INTERPTFR interpolate time and frequency bins of the time-frequency spectrum data.
%
% Use as
%   cfg           = [];
%   cfg.parameter = 'fourierspctrm';
%   cfg.method    = 'spline';
%   data = ft_interpTFR(cfg,data)
%
% The input freq structure should be a a time-frequency representation of
% power or coherence that was computed using the FT_FREQANALYSIS function.
% 'keeptrials' data is also valid, but calculation time may be long.
% 
% If data contains a lot of NaN values, specifying exptrapval and baseline
% is recommended or 'spline' method. (Example: cfg.extrapval = 0; cfg.baseline = 'yes';)
% 
% The baseline correction is applied before the interpolation. If you want
% to compute the inpterpolation on raw data, please use 'cfg.baseline'
% option in 'ft_singleplotTFR',  'ft_topoplotTFR' and 'ft_multiplotTFR'.
%
% The configuration can have the following parameters:
%   cfg.parameter      = field to be plotted on z-axis, e.g. 'powspcrtrm' (default depends on data.dimord)
%   cfg.method         = 'nearest', 'linear', 'cubic', or 'spline' (default = 'linear')
%                        different methods of calculating the interpolation, see INTERP2
%   cfg.extrapval      = scalar or 'none' (default = 'none')
%                        function value outside domain of the input data, see INTERP2
%                        NaN values in the input data are replaced by this value
%   cfg.keepnan        = 'yes' or 'no' (default = 'yes')
%                        if 'no', fill interpolated data into outside domain of input data
%   cfg.tintvl         = number, new interval of time dimension (default = 'auto')
%                        if 'auto', set the minimum interval in 'data' 
%   cfg.fintvl         = number, new interval of frequency dimension (default = 'auto')
%                        if 'auto', set the minimum interval in 'data' 
%   cfg.toilim         = [begin end], time band of interest (default = 'all')
%                        new toi is [begin cfg.tintvl end]
%   cfg.foilim         = [begin end], frequency band of interest (default = 'all')
%                        new foi is [begin cfg.fintvl end]
% 
% The following options are available:
%   cfg.baseline       = 'yes','no' or [time1 time2] (default = 'no'), see FT_FREQBASELINE
%   cfg.baselinetype   = 'absolute', 'relative', 'relchange', 'vssum' or 'db' (default = 'absolute')
% 
% Undocumented options: 
%   cfg.toi
%   cfg.foi
%       both two options are refered instead of data.time and data.freq when cfg.*intvl or cfg.*oilim is 'auto'
% 
%   See also FT_FREQANALYSIS, FT_FREQINTERPOLATE, FT_FREQBASELINE, FT_FOI_ADJUST, INTERP2.

%   using: fieldTrip, ft_private(getdimord), cellstrfind, nearlyeq

%   20170205 Yuasa
%   20170221 Yuasa: add baseline correction
%   20170607 Yuasa: enable cfg.toi, cfg.foi
%   20170921 Yuasa: separate computation method for complex matrix & real matrix
%   20190322 Yuasa: minor update

ft_defaults;
ft_private('silent');

%-- check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'freq');

%-- Set the defaults:
cfg.parameter      = ft_getopt(cfg, 'parameter', 'powspctrm');
cfg.method         = ft_getopt(cfg, 'method', 'linear');
cfg.extrapval      = ft_getopt(cfg, 'extrapval', 'none');
cfg.keepnan        = ft_getopt(cfg, 'keepnan', 'yes');
cfg.tintvl         = ft_getopt(cfg, 'tintvl', 'auto');
cfg.fintvl         = ft_getopt(cfg, 'fintvl', 'auto');
cfg.toilim         = ft_getopt(cfg, 'toilim', 'all');
cfg.foilim         = ft_getopt(cfg, 'foilim', 'all');
cfg.baseline       = ft_getopt(cfg, 'baseline',     'no');
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.maskparameter  = ft_getopt(cfg, 'maskparameter', []);

if isempty(cfg.tintvl),   cfg.tintvl  = 'auto';     end
if isempty(cfg.fintvl),   cfg.tintvl  = 'auto';     end
if isempty(cfg.toilim) || strcmp(cfg.toilim,'auto'),    cfg.toilim  = 'all';     end
if isempty(cfg.foilim) || strcmp(cfg.foilim,'auto'),    cfg.foilim  = 'all';     end

%-- check spacing uniformity
if strcmp(cfg.method,'cubic')
  timespace   = abs(diff(data.time));
  freqspace   = abs(diff(data.freq));
  if numel(unique(timespace))>1 || numel(unique(freqspace))>1
    if (maxmindiff(timespace) > min(timespace) * 10^-3) || (maxmindiff(freqspace) > min(freqspace) * 10^-3)       % if difference of difference is larger than 0.1% of difference
      warning('%s\n%s','The ''cubic'' method requires the grid to have a uniform spacing.','Switching the method from ''cubic'' to ''spline'' because this condition is not met.');
      cfg.method    = 'spline';
    end
  end
end

%-- check dimord
if isfield(data, 'dimord'),
 dimord = data.dimord;
else
 dimord = getdimord(data, cfg.parameter);
end
dimtok = tokenize(dimord, '_');
dimtime = cellstrfind(dimtok, 'time');
dimfreq = cellstrfind(dimtok, 'freq');

if ~isfield(data, cfg.parameter)
  error('data has no field ''%s''', cfg.parameter);
end
if ndims(data.(cfg.parameter)) ~= length(dimtok)
  error('fail to specify the structure of ''%s''', cfg.parameter);
elseif isempty(dimtime) || isempty(dimfreq)
  error('input data needs a time dimension and frequency dimension');
end

dimother = 1:length(dimtok); dimother([dimfreq dimtime]) = [];
dimperm  = 1:length(dimtok);
dimperm(2,[dimother dimfreq dimtime]) = 1:length(dimtok);
dimperm  = sortrows(dimperm',2)';

%-- apply baseline correction
if ~strcmp(cfg.baseline, 'no')
  %-- keep mask-parameter if it is set
  if ~isempty(cfg.maskparameter)
    tempmask = data.(cfg.maskparameter);
  end
  data = ft_freqbaseline(cfg, data);
  %-- put mask-parameter back if it is set
  if ~isempty(cfg.maskparameter)
    data.(cfg.maskparameter) = tempmask;
  end
end

%-- reshape spectrum into (:,freq,time)
% dimtokperm  = dimtok(dimperm(1,:));
spctrm      = permute(data.(cfg.parameter),dimperm(1,:));
nlength     = size(spctrm);
spctrm      = reshape(spctrm,[],nlength(end-1),nlength(end));

%-- set 'auto' val
if ischar(cfg.tintvl) && strcmp(cfg.tintvl, 'auto')
  if isfield(cfg,'toi'),    cfg.tintvl   = min(diff(cfg.toi));
  else                      cfg.tintvl   = min(diff(data.time));
  end
end
if ischar(cfg.fintvl) && strcmp(cfg.fintvl, 'auto')
  if isfield(cfg,'foi'),    cfg.fintvl   = min(diff(cfg.foi));
  else                      cfg.fintvl   = min(diff(data.freq));
  end
end
if ischar(cfg.toilim) && strcmp(cfg.toilim, 'all')
  if isfield(cfg,'toi'),    cfg.toilim   = [min(cfg.toi) max(cfg.toi)];
  else                      cfg.toilim   = [min(data.time) max(data.time)];
  end
end
if ischar(cfg.foilim) && strcmp(cfg.foilim, 'all')
  if isfield(cfg,'foi'),    cfg.foilim   = [min(cfg.foi) max(cfg.foi)];
  else                      cfg.foilim   = [min(data.freq) max(data.freq)];
  end
end

if ~isnumeric(cfg.tintvl) || ~isnumeric(cfg.fintvl)
    error('interval parameter is invalid');
end
if ~isnumeric(cfg.toilim) || ~isnumeric(cfg.foilim) || numel(cfg.toilim)<2 || numel(cfg.foilim)<2
    error('toi or foi is invalid');
end

%%%% main %%%%
tic;
%-- adjust toilim & foilim (check whether distance to the nearest value is less than 1% of difference)
nearidx     = nearlyeq(data.time, cfg.toilim);
if cfg.toilim(1) < data.time(nearidx(1)) && ...
        abs(diff([cfg.toilim(1), data.time(nearidx(1))])) < min([abs(diff(data.time)), cfg.tintvl(1)])*0.01,
    data.time(nearidx(1)) = cfg.toilim(1);
end
if cfg.toilim(end) > data.time(nearidx(end)) && ...
        abs(diff([cfg.toilim(end), data.time(nearidx(end))])) < min([abs(diff(data.time)), cfg.tintvl(end)])*0.01,
    data.time(nearidx(end)) = cfg.toilim(end);
end
nearidx     = nearlyeq(data.freq, cfg.foilim);
if cfg.foilim(1) < data.freq(nearidx(1)) && ...
        abs(diff([cfg.foilim(1), data.freq(nearidx(1))])) < min([abs(diff(data.freq)), cfg.tintvl(1)])*0.01,
    data.freq(nearidx(1)) = cfg.foilim(1);
end
if cfg.foilim(end) > data.freq(nearidx(end)) && ...
        abs(diff([cfg.foilim(end), data.freq(nearidx(end))])) < min([abs(diff(data.freq)), cfg.tintvl(end)])*0.01,
    data.freq(nearidx(end)) = cfg.foilim(end);
end

%-- get original grid
[gridtime, gridfreq]=meshgrid(data.time,data.freq);

%-- make new grid
newtime = cfg.toilim(1):cfg.tintvl(1):cfg.toilim(end);
newfreq = cfg.foilim(1):cfg.fintvl(1):cfg.foilim(end);

[gridnewtime, gridnewfreq]=meshgrid(newtime, newfreq);

%-- interpolation loop
warn_state = warning('off','MATLAB:interp2:NaNstrip');      % skip NaNstrip warning
try
newspctrm   = zeros(size(spctrm,1),length(newfreq), length(newtime));
if strcmp(cfg.keepnan,'yes')
    ispctrm 	= squeeze(mean(spctrm(:,:,:),1));
    nanspctrm   = isnan(interp2(gridtime,gridfreq, ispctrm, gridnewtime,gridnewfreq, 'nearest'));
end
for iloop = 1:size(spctrm,1)

    %-- get spectrum
    ispctrm 	= squeeze(spctrm(iloop,:,:));

    %-- fill 0 at NaN
    if ~strcmp(cfg.extrapval,'none')
      ispctrm(isnan(ispctrm)) = cfg.extrapval;
    end

    %-- interpolation
    if isreal(ispctrm)      % for real matrix
     if strcmp(cfg.extrapval,'none')
    newspctrm(iloop,:,:)   = interp2(gridtime,gridfreq, ispctrm, gridnewtime,gridnewfreq, cfg.method);
     else
    newspctrm(iloop,:,:)   = interp2(gridtime,gridfreq, ispctrm, gridnewtime,gridnewfreq, cfg.method, cfg.extrapval);
     end
    else                    % for complex matrix
      %-- power: interpolate in power-space & delete negative value
     if strcmp(cfg.extrapval,'none')
    newspctrmpow   = interp2(gridtime,gridfreq, abs(ispctrm), gridnewtime,gridnewfreq, cfg.method);
     else
    newspctrmpow   = interp2(gridtime,gridfreq, abs(ispctrm), gridnewtime,gridnewfreq, cfg.method, cfg.extrapval);
     end
    newspctrmpow(newspctrmpow<0) = 0;
      %-- phase: interpolate in complex-space & convert into phase
     if strcmp(cfg.extrapval,'none')
    newspctrmphs   = interp2(gridtime,gridfreq, ispctrm, gridnewtime,gridnewfreq, cfg.method);
     else
    newspctrmphs   = interp2(gridtime,gridfreq, ispctrm, gridnewtime,gridnewfreq, cfg.method, cfg.extrapval);
     end
    newspctrmphs   = newspctrmphs ./ abs(newspctrmphs);
    newspctrmphs(isnan(newspctrmphs)) = 0;
      %-- combine power and phase
    newspctrm(iloop,:,:)                = newspctrmpow .* newspctrmphs;
    end
    
    %-- re-fill NaN
    if strcmp(cfg.keepnan,'yes')
      newspctrm(iloop,nanspctrm)  = nan;
    end
end
catch ME
 ft_private('remove','silent');
 warning(warn_state);
 if strcmp(ME.identifier,'MATLAB:interp2:NotEnoughPointsNanStrip')
    itpEr = [];
    itpEr.identifier = ME.identifier;
    itpEr.message = sprintf('%s\nPlease modify input data or specify cfg.extrapval.',ME.message);
    itpEr.stack   = ME.stack(end);
    error(itpEr);
 else
    rethrow(ME);
 end
end

%-- return reshape & re-permute
newspctrm   = reshape(newspctrm, [nlength(1:end-2), length(newfreq), length(newtime)]);
dimperm     = sortrows(dimperm',1)';
newspctrm   = permute(newspctrm,dimperm(2,:));

%-- output
data.time  = newtime;
data.freq  = newfreq;
data.(cfg.parameter) = newspctrm;

fprintf('the call to  ''ft_interpTFR'' took %.0f seconds to interpolate the spectrum data\n',toc);

ft_private('remove','silent');
warning(warn_state);

end

function dd = maxmindiff(d)
% get difference of maximum and minimum

dd = diff([min(d(:)), max(d(:))]);

end

%%% check script
%{
figure
surf(gridnewtime,gridnewfreq, abs(newspctrm))
view([0 90]); grid off; clim([0 2.5]); xlim([-0.5 2]); ylim([5 80]); colormap('jet')
set(get(gca,'Children'),'LineStyle','none')
set(gcf,'MenuBar','none','ToolBar','none')

figure
surf(gridnewtime,gridnewfreq, angle(newspctrm))
view([0 90]); grid off; clim([-pi pi]); xlim([-0.5 2]); ylim([5 80]); colormap('hsv')
set(get(gca,'Children'),'LineStyle','none')
set(gcf,'MenuBar','none','ToolBar','none')
%}