function [data] = ft_regressout(cfg,data)
% FT_REGRESSOUT regresses out event related responses from single trials time seriese. 
% When cfg.maxlag is specified, time lags of event related responses across trials 
% are estimated based on cross-correlation.
% 
% Use as
%   [data_out] = ft_regressout(cfg, data)
% 
% The configuration can have the following parameters:
%   cfg.demean         = 'no' or 'yes', whether to apply baseline correction (default = 'no')
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial (default = 'all')
%   cfg.groups         = vector 1 x trials, consider trials having same index as the same condition,
%                        the default considers all trials as the same condition (default = [])
%   cfg.maxlag         = maximum time lag in seconds (default = 0)
%   cfg.negregress     = no' or 'yes', whether to allow negative regression during time lag estimation
%                       (default = 'no')
%   cfg.refwindow      = [begin end] in seconds, evaluate time lag in the reference latencies, 
%                        the default is the complete trial (default = 'all')
%     * too short reference window may cause regression to a local maxima of noise signals .
%       shold include pre adn post baseline period (eg.[-0.2 0.6]).
%
% You can check time series of regressor as:
%   regressor = cellfun(@(x,y) x-y, data.trial, data_out.trial, 'UniformOutput', false);

% 20190807 Yuasa
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig
end

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set private functions
ft_private('set','silent')
try

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% set the defaults
cfg.trials         = ft_getopt(cfg, 'trials',        'all');
cfg.channel        = ft_getopt(cfg, 'channel',       'all');
cfg.demean         = ft_getopt(cfg, 'demean',        'no');
cfg.baselinewindow = ft_getopt(cfg, 'baselinewindow','all');
cfg.groups         = ft_getopt(cfg, 'groups',         []);
cfg.maxlag         = ft_getopt(cfg, 'maxlag',         0);
cfg.negregress     = ft_getopt(cfg, 'negregress',    'no');
cfg.refwindow      = ft_getopt(cfg, 'refwindow',     'all');

isdemean        = strcmp(cfg.demean, 'yes');
isnegregress    = strcmp(cfg.negregress, 'yes');

% select trials of interest
tmpcfg = [];
tmpcfg.trials = cfg.trials;
tmpcfg.channel = cfg.channel;
data = ft_selectdata(tmpcfg, data);
% restore the provenance information

[cfg, data] = rollback_provenance(cfg, data);

% check time length
begsamplatency = cellfun(@min,data.time);
endsamplatency = cellfun(@max,data.time);
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];
assert(all(minperlength==maxperlength),'data has variable trial lengths, you specified not to accept that');

% check group list
if isempty(cfg.groups)
    cfg.groups  = ones(1,length(data.trial));
end
assert(length(cfg.groups)==length(data.trial),'cfg.groups must have the same number of trials');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main
fsample     = data.fsample;
time        = data.time{1};
ntrial      = length(data.trial);
nchan       = length(data.label);   % number of channels
nsamples    = length(time);

maxlag      = cfg.maxlag .* fsample;

% baseline correction
if isdemean
    if ischar(cfg.baselinewindow) && strcmp(cfg.baselinewindow, 'all')
      % the begin and endsample of the baseline period correspond to the complete data minus padding
      blcsample = 1:nsamples;
    else
      % determine the begin and endsample of the baseline period and baseline correct for it
      blcsample = nearest(time, cfg.baselinewindow(1)):nearest(time, cfg.baselinewindow(2));
    end
    data.trial       = cellfun(@(dat) bsxfun(@minus,dat,mean(dat(:,blcsample),2)),data.trial,'UniformOutput',false);
end

% reference latencies for time lag estimation
if ischar(cfg.refwindow) && strcmp(cfg.refwindow, 'all')
    lagrefsample = 1:nsamples;
else
    lagrefsample = nearest(time, cfg.refwindow(1)):nearest(time, cfg.refwindow(2));
end


% regress erp out
dat     = cat(3,data.trial{:});
dat_out = nan(nchan,nsamples,ntrial);
condlist = reshape(unique(cfg.groups),1,[]);
param_reg = zeros(ntrial,nchan);
param_lag = zeros(ntrial,nchan);
for ich=1:nchan
    fprintf('regressing erp from channel %d\n', ich);
    for icond=condlist
        % ERP model
        curcond = cfg.groups==icond;
        av_erp  = reshape(squeeze(mean(dat(ich,:,curcond),3)),nsamples,1);
        if all(diff(av_erp)==0)
            % skip regress if regressor is flat
            dat_out(ich,:,curcond) = reg_R(ich,:,curcond);
        else
            for iepoch=find(curcond)
                x  = reshape(squeeze(dat(ich,:,iepoch)),nsamples,1);
                % estimate time lag
                [cr,tlags] = xcorr(x(lagrefsample),av_erp(lagrefsample),maxlag);
                if all(diff(cr)==0)
                    tlags = 0; dt = 1;  % avoid empty output for flat cr
                elseif isnegregress
                    [~,dt] = findpeaks(abs([0 cr' 0]).^2,'Npeaks',1,'SortStr','descend');
                    dt = dt - 1;        % add 0 at both side of cr to avoid epmty output
                else
                    [~,dt] = findpeaks([min(cr) cr' min(cr)],'Npeaks',1,'SortStr','descend');
                    dt = dt - 1;        % add 0 at both side of cr to avoid epmty output
                end
                % regress ERP out
                reg_erp = zeros(nsamples,1);
                reg_erp(max(1,tlags(dt)+1):min(nsamples,nsamples+tlags(dt))) = ...
                    av_erp(max(1,-tlags(dt)+1):min(nsamples,nsamples-tlags(dt)));
                [B,~,reg_R] = regress(x,reg_erp);
                dat_out(ich,:,iepoch) = reg_R;
                param_reg(iepoch,ich) = B;
                param_lag(iepoch,ich) = tlags(dt) ./ fsample;
            end
        end
    end
end
data_out = squeeze(mat2cell(dat_out,nchan,nsamples,ones(ntrial,1)))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data

data.trial    = data_out;
data.regressB   = param_reg;
data.regressLag = param_lag;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

catch ME
    ft_private('reset','silent')
    rethrow(ME)
end % try
