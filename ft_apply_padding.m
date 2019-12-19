function    [data] = ft_apply_padding(cfg, data)

% FT_APPLY_PADDING  performs padding on the data, i.e. adds or removes samples
% to or from the fieldTrip data.
%
% Use as
%   [data] = ft_apply_padding(cfg, data)
% 
% The configuration can have the following parameters:
%   cfg.padding      = length (in seconds) to which the trials are padded for filtering (default = 0)
%                      data is padded until the time length reaches to this value
%   cfg.padtype      = 'zero', 'mean', 'localmean', 'edge', 'mirror', 'nan' or 'remove' (default = 'zero')
%                      or {pre_padtype, post_padtype} is valid
%   cfg.paddir       = direction of padding, 'left'/'right'/'both' (default = 'both')
% 
% See also FT_PREPROC_PADDING, FT_PREPROCESSING

% using: fieldtrip

% 20170220 Yuasa: pickup from ft_preprocessing and modified

narginchk(2,2);

ft_defaults;

% set all the defaults
cfg.padding         = ft_getopt(cfg, 'padding',      0);
cfg.paddir          = ft_getopt(cfg, 'paddir',  'both');
cfg.padtype         = ft_getopt(cfg, 'padtype', 'zero');

if ~iscell(cfg.padtype),    cfg.padtype = {cfg.padtype};    end

%-- check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'mvar'});
%-- determine sampling frequency
fsample = data.fsample;


%-- set padding
if cfg.padding>0
    padding = round(cfg.padding * fsample);
    %-- update the configuration (in seconds) for external reference
    cfg.padding = padding / fsample;
else
    %-- no padding was requested
    padding = 0;
end

%-- loop for padding
ntrl = length(data.trial);
for itrl=1:ntrl
    nsamples = numel(data.time{itrl});
    
    % pad data by mirroring
    if nsamples>padding
      % the trial is already longer than the total length requested
      begpadding = 0;
      endpadding = 0;
    else
      switch cfg.paddir
        case 'both'
          % begpadding+nsamples+endpadding = total length of data
          begpadding = ceil((padding-nsamples)/2);
          endpadding = floor((padding-nsamples)/2);
        case 'left'
          begpadding = padding-nsamples;
          endpadding = 0;
        case 'right'
          begpadding = 0;
          endpadding = padding-nsamples;
        otherwise
          error('unsupported requested direction of padding');
      end
    end
    begtime = data.time{itrl}(1)   - begpadding / fsample;
    endtime = data.time{itrl}(end) + endpadding / fsample;
    
    if length(cfg.padtype)>1
      data.trial{itrl} = ft_preproc_padding(data.trial{itrl}, cfg.padtype{1}, begpadding, 0);
      data.trial{itrl} = ft_preproc_padding(data.trial{itrl}, cfg.padtype{2}, 0, endpadding);
    else
      data.trial{itrl} = ft_preproc_padding(data.trial{itrl}, cfg.padtype{1}, begpadding, endpadding);
    end
      data.time{itrl}  = [begtime:1/fsample:data.time{itrl}(1) data.time{itrl}(2:end-1) data.time{itrl}(end):1/fsample:endtime];
    
end % for all trials

end