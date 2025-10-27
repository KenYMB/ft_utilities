function [phase] = ft_phaseanalysis(cfg, freq)

% FT_PHASEANALYSIS computes inter-trial phase locking.
%
% Use as
%   phase = ft_connectivityanalysis(cfg, freq)
% 
% The configuration structure has to contain
%   cfg.method  =  string, can be
%     'itpc',        inter-trial phase coherence (default)
%                    inter-trial phase-locking factor
%     'itlc',        inter-trial linear coherence
%                    inter-trial coherence
%
% Additional configuration options are
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.frequency     = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.latency       = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')

% Using: ft_private

% 20190318 Yuasa
% 20251024 Yuasa: add try for ft_preamble to avoid compatibility errors

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
try
ft_preamble init
ft_preamble debug
ft_preamble loadvar freq
ft_preamble provenance freq
ft_preamble trackconfig
end

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.method     = ft_getopt(cfg, 'method',     'itpl');
cfg.trials     = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel    = ft_getopt(cfg, 'channel',    'all');
cfg.frequency  = ft_getopt(cfg, 'frequency',  'all');
cfg.latency    = ft_getopt(cfg, 'latency',    'all');

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', {'freq', 'freqmvar'}, 'feedback', 'yes');
inparam = 'fourierspctrm';

% determine some specific details of the input data
hasphs   = isfield(freq, inparam);
hasrpt   = ~isempty(strfind(freq.dimord, 'rpt')) || ~isempty(strfind(freq.dimord, 'subj'));
hastim   = ~isempty(strfind(freq.dimord, 'time'));

% check sensibility of configuration
if ~hasphs, error('inter-trial phase analysis requires input data with ''fourierspctrm''');                                                end
if ~hasrpt, error('inter-trial phase analysis requires input data with repeated observations');                                            end

% select data of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'latency', 'frequency'});
freq = ft_selectdata(tmpcfg, freq);
% restore the provenance information
ft_private('silent')
[cfg, freq] = rollback_provenance(cfg, freq);
ft_private('reset','silent')

%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%
F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);            % number of trials

switch cfg.method
  case {'itpl','itplf','plf','itpc'}
      outparam      = 'itpcspctrm';
      datout        = F./abs(F);          % divide by amplitude
      datout        = nanmean(datout,1);  % sum angles and normalize
      datout        = abs(datout);        % take the absolute value
      datout        = squeeze(datout);    % remove the first singleton dimension
  case {'itcoh','coh','itlc','itc'}
      outparam      = 'itlcspctrm';
      datout        = nansum(F) ./ (sqrt(N*nansum(abs(F).^2)));
      datout        = abs(datout);        % take the absolute value and normalize, i.e. ignore phase
      datout        = squeeze(datout);    % remove the first singleton dimension
  otherwise
    error('unknown method % s', cfg.method);
end

% create the output structure
phase = [];
phase.label = freq.label;
phase.(outparam) = datout;
if isfield(freq, 'dimord'),
    dimtok = tokenize(freq.dimord, '_');
else
    ft_private('silent');
    dimtok = tokenize(getdimord(freq, inparam), '_');
    ft_private('reset','silent')
end
dimord = '';
for k = 1:numel(dimtok)
    if isempty(strfind(dimtok{k}, 'rpt'))
        dimord = [dimord, '_', dimtok{k}];
    end
end
phase.dimord = dimord(2:end);
if any(strcmp(dimtok, 'time')), phase.time = freq.time; end
if any(strcmp(dimtok, 'freq')), phase.freq = freq.freq; end
if isfield(freq, 'grad'), phase.grad = freq.grad; end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   freq
ft_postamble provenance phase
ft_postamble history    phase
ft_postamble savevar    phase
