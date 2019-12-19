function ft_neurofif_write_average(cfg,varargin)
%
% FT_NEUROFIF_WRITE_AVERAGE(cfg,timelock1,timelock2,...) converts fieldtrip
% timelocked average data into neuromag evoked data, and saves fiff file.
% 
% The configuration should be according to
% 
%   cfg.output     = output file name
%   cfg.rawdata    = file name of original fiff file (header information of output refers this file)
%                    'auto': load data.cfg.((previous)).dataset (default)
%   cfg.badch      = Nx1 cell-array with selection of bad channels (default = ''),
%                    'mag'/'grad'/'vert': select magnetmeter/planar/vertical sensor as bad
%                    '*' is interpreted as a wildcard
%   cfg.datalabel  = Nx1 cell-array with comment for each timelocked average data
%                    if empty, input order is applied (default)
%   cfg.fsample    = sampling frequency of the averaged data [Hz]
%                    if empty, inverse the periodic time of input data
%   cfg.unitcorr   = yes' or 'no', to correct magnetometer unit into 'T'
%                    and gradient sensor unit into 'T/m' (default = 'no')
%   *** Caution ***
%   If you load MEG data with ft_preprocessing, the sensor units for
%   planars are considered as T/cm while the intact units for planar
%   sensors are T/m. cfg.unitcorr is only usefull if you corrected this
%   contradiction manually.
% 
% Covariance file is also output if timelocked average data include covariance matrix.
% The configuration should be according to (See also ft_mnefif_write_covariance)
% 
%   cfg.cov.output  = output file name for covariance data
%                     if empty, saving covariance file is skipped
%   cfg.cov.merge   = 'yes' or 'no', merge covariance across all input data
%                     correct the degrees of freedom before export (default = 'yes')
%   cfg.cov.dofcorr = 'yes' or 'no', to correct the computational method of
%                     defrees of freedom of covariance (default = 'yes')
% ( cfg.cov.diag    = yet not support)
% 
% Undocumented option
% 
%   cfg.overwrite  = 'yes' or 'no', (default = 'no')
% 
% 
% Exapmle1: 
%   cfg            = [];
%   cfg.output     = 'MEG_ave.fif';
%   cfg.badch      = 'vert';
%   cfg.datalabel  = 'Cond_1';
%   cfg.unitcorr   = 'yes';                 % reset manual correction
%   ft_neurofif_write_average(cfg,avgMEG);
% 
% Exapmle2: 
%   cfg            = [];
%   cfg.output     = 'MEG_ave.fif';
%   cfg.badch      = 'vert';
%   cfg.datalabel  = {'Cond_1','Cond_2'};
%   cfg.unitcorr   = 'yes';                 % reset manual correction
%   ft_neurofif_write_average(cfg,avgMEG1,avgMEG2);
% 
% See also FIFF_WRITE_EVOKED, FT_NEUROFIF_WRITE_RAW, FT_MNEFIF_WRITE_COVARIANCE

% using: fieldtrip(, ft_hastoolbox, mne, match_str), cellstrfind, int2ordinal

% 20170821 Yuasa
% 20170829 Yuasa: separate covariance output parts
% 20190508 Yuasa: enable '-' option in cfg.badch


curpath = path;
ft_defaults;
ft_hastoolbox('mne',1,1);

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
me = 'original:ft_neurofif_write_average';
if nargin < 2
    path(curpath);
    error(me,'Incorrect number of arguments');
end

try
%--  prepare cfg
cfg.rawdata     = ft_getopt(cfg,'rawdata','auto');
cfg.output      = ft_getopt(cfg,'output','dummy_ave.fif');
cfg.badch       = ft_getopt(cfg,'badch');
cfg.datalabel   = ft_getopt(cfg,'datalabel');
cfg.fsample     = ft_getopt(cfg,'fsample');
cfg.unitcorr    = ft_getopt(cfg,'unitcorr','no');
cfg.cov         = ft_getopt(cfg,'cov');
cfg.cov.output  = ft_getopt(cfg.cov,'output','');
cfg.cov.merge   = ft_getopt(cfg.cov,'merge','yes');
cfg.cov.dofcorr = ft_getopt(cfg.cov,'dofcorr','yes');
cfg.cov.diag    = ft_getopt(cfg.cov,'diag',false);
cfg.overwrite   = ft_getopt(cfg,'overwrite','no');
%--- bad channels
if ~isempty(cfg.badch)                                                  	% interpret cfg.badch
    if  ~iscell(cfg.badch)
        cfg.badch = {cfg.badch};
    else
        cfg.badch = reshape(cfg.badch,1,[]);
    end
    genidx    = match_str(cfg.badch,{'all','mag','grad','vert'});
    for ich = 1:length(genidx)
        switch cfg.badch{genidx(ich)}
            case 'all'
                cfg.badch = [cfg.badch {'MEG*'}];
            case 'mag'
                cfg.badch = [cfg.badch {'MEG*1' 'MEG*4'}];
            case 'grad'
                cfg.badch = [cfg.badch {'MEG*2' 'MEG*3' 'MEG*5' 'MEG*6'}];
            case 'vert'
                cfg.badch = [cfg.badch {'MEG*4' 'MEG*5' 'MEG*6'}];
        end
    end
    cfg.badch(genidx) = [];
end
%--- label for each data
if ~isempty(cfg.datalabel) && ~iscell(cfg.datalabel)
    cfg.datalabel = {cfg.datalabel};
end
if length(cfg.datalabel) < (nargin - 1)
    labelparts    = cfg.datalabel;
    cfg.datalabel = cellstr(num2str([1:(nargin - 1)]','%-d'));
    cfg.datalabel(1:length(labelparts)) = labelparts;
end
%--- unit correction
isucorr  = strcmpi(cfg.unitcorr,'yes');
%--- covdat flag
iscov    = ~isempty(cfg.cov.output);                                        % update in loop
%--- overwrite flag
isforce  = strcmpi(cfg.overwrite,'yes');

%-- check raw-data
TMPdata = varargin{1};
datasetcand = ft_findcfg(TMPdata.cfg, 'datafile');
if strcmp(cfg.rawdata,'auto')
    if ~isempty(datasetcand);
        cfg.rawdata  = datasetcand;
    else
        error(me,'Failed to detect raw file. Please specify cfg.rawdata.');
    end
end

%-- get information of raw-data
data = [];
if exist(cfg.rawdata, 'file')
    %--- copy field
    try     data.info = fiff_read_meas_info(cfg.rawdata);
    catch,  error(me,'File %s is invalid',cfg.rawdata);
    end
else
    error(me,'Cannot open file %s',cfg.rawdata);
end
%--- check field
if ~isfield(data.info,'meas_id'),       data.info.meas_id    = struct([]);  end
if ~isfield(data.info,'meas_date'),     data.info.meas_date  = [];          end
if ~isfield(data.info,'projs'),         data.info.projs      = struct([]); 	end
if ~isfield(data.info,'comps'),         data.info.comps      = struct([]);  end
%--- following info is copied directly from raw-data in 'fiff_write_evoked'
if ~isfield(data.info,'dev_head_t'),    data.info.dev_head_t = struct([]);  end     % transform matrix
if ~isfield(data.info,'ctf_head_t'),    data.info.ctf_head_t = struct([]);  end
if ~isfield(data.info,'dig'),           data.info.dig        = struct([]);  end     % digitizer poitns

%-- modify channel info
if isfield(TMPdata,'grad')
    labelmatch = match_str(data.info.ch_names, TMPdata.grad.label);
else
    labelmatch = match_str(data.info.ch_names, TMPdata.label);
end
data.info.chs       = data.info.chs(labelmatch);
data.info.ch_names  = data.info.ch_names(labelmatch);
data.info.nchan     = length(data.info.chs);

%-- input preprocessing parameter
if isempty(cfg.fsample)
  data.info.sfreq     = nanmin(round(1./diff(TMPdata.time)));
else
  data.info.sfreq     = cfg.fsample;
end
data.info.highpass  = nanmax([data.info.highpass, nanmin(ft_findcfg(TMPdata.cfg,'bpfreq')), ft_findcfg(TMPdata.cfg,'hpfreq')]);
data.info.lowpass   = nanmin([data.info.lowpass, nanmax(ft_findcfg(TMPdata.cfg,'bpfreq')), ft_findcfg(TMPdata.cfg,'lpfreq')]);

%-- time-locked data
nanchs = [];
data.evoked(1:(nargin-1)) = struct('aspect_kind',FIFF.FIFFV_ASPECT_AVERAGE,'is_smsh',0,'nave',1,'first',0,'last',0,'comment','','times',[],'epochs',[]);
for idats = 1:(nargin-1)
    %--- check data
    TMPdata = varargin{idats};
    assert(isfield(TMPdata,'avg'),me,'The %s data does not include averaged data', int2ordinal(idats));
    if ~strcmp(datasetcand, ft_findcfg(TMPdata.cfg, 'datafile'))
        warning(me,'input data comes from different raw-datafiles');
        cfg.info.rawdata   = '';              % remove rawdata information
    end
        
    %-- check sample
    sampletime = round(diff(TMPdata.time).*data.info.sfreq)./data.info.sfreq;
    if (min(sampletime(:)) ~= max(sampletime(:))) || (round(1/sampletime(1)) ~= data.info.sfreq)
        error(me, 'The sampling time points in %s data is inconsistent with the sampling frequency', int2ordinal(idats));
    end
    
    %-- unit correction
    if isucorr
      if isfield(TMPdata,'grad')
         %--- convert fT into T
         labelmatch     = match_str(TMPdata.grad.label, TMPdata.label);
         corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),{'fT','fT*'});
         TMPdata.avg(corrchan,:)    = TMPdata.avg(corrchan,:) * 1e-15;
         %--- convert mm/cm/dm into m
         corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),'*/mm');
         TMPdata.avg(corrchan,:)    = TMPdata.avg(corrchan,:) / 1e-3;
         corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),'*/cm');
         TMPdata.avg(corrchan,:)    = TMPdata.avg(corrchan,:) / 1e-2;
         corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),'*/dm');
         TMPdata.avg(corrchan,:)    = TMPdata.avg(corrchan,:) / 1e-1;
      else
         warning(me, 'Cannot specify current units, ''grad'' field does not exist');
      end
    end
    
    %-- first and last samples
    data.evoked(idats).first     = round(TMPdata.time(1) * data.info.sfreq);
    data.evoked(idats).last      = round(TMPdata.time(end) * data.info.sfreq);
    
    %-- bad channel update
    tmpbads = min(isnan(TMPdata.avg),[],2);
    TMPdata.avg(tmpbads,:) = 0;
    nanchs  = [nanchs, reshape(TMPdata.label(tmpbads),1,[])];
    
    %-- timelock data (check channel info & bad channel)
    [labelmatch, labelidx] = match_str(data.info.ch_names, TMPdata.label);
    data.evoked(idats).epochs                = zeros(length(data.info.ch_names), size(TMPdata.avg,2));
    data.evoked(idats).epochs(labelmatch,:)  = TMPdata.avg(labelidx,:);
    
    %-- number of trials
    if isfield(TMPdata,'dof')
        data.evoked(idats).nave  = nanmin(TMPdata.dof(:));
    else
        warning(me, 'Failed to find the number of averaged trials in %s data', int2ordinal(idats));
    end
    
    %-- comment
    data.evoked(idats).comment = cfg.datalabel{idats};
    
    %-- covariance
    if iscov && ~isfield(TMPdata,'cov')
        warning(me,'The %s data does not include covariance matrix\nSkip saving covariance data', int2ordinal(idats));
        iscov = false;
    end
    
end

%-- input preprocessing parameter
negidx           = ~cellfun(@isempty,regexp(cfg.badch,'^-','once'));
negchs   = cellstrfind(data.info.ch_names, regexprep(cfg.badch(negidx),'^-',''));
    negchs       = idx2logical(data.info.ch_names,negchs);
badchs   = cellstrfind(data.info.ch_names, cfg.badch(~negidx));
    badchs       = idx2logical(data.info.ch_names,badchs);
badchs   = badchs & ~negchs;                                                % exclude '-*' chs from bad-channels
cfg.badch        = data.info.ch_names(badchs);                            	% register bad ch
nanchs   = cellstrfind(data.info.ch_names, nanchs);
    nanchs       = idx2logical(data.info.ch_names,nanchs);
data.info.bads   = data.info.ch_names(badchs | nanchs);

%%%
%-- prepare for output
[savedir,savename,saveext] = fileparts(cfg.output);         % check extention
if ~strcmp(saveext,'.fif') && ~strcmp(saveext,'.fiff')
    cfg.output = [cfg.output '.fif'];
end
if ~isforce && logical(exist(cfg.output,'file'))            % check overwriting
    curwarn = warning('off', 'backtrace');
    warning('The output file for average already exist! Now copying...');
    warning(curwarn);
    copyfile(cfg.output, [cfg.output '.bk']);
end
if iscov
  [savedir,savename,saveext] = fileparts(cfg.cov.output);	% check extention
  if ~strcmp(saveext,'.fif') && ~strcmp(saveext,'.fiff')
    cfg.cov.output = [cfg.cov.output '.fif'];
  end
  [savedir,savename,saveext]    = fileparts(cfg.output);    % when the names are same
  [savedirc,savenamec,saveextc] = fileparts(cfg.cov.output);
  if strcmpi(savename, savenamec)
    curwarn = warning('off', 'backtrace');
    warning('The output file for average and for covariance are same! Changing the name...');
    warning(curwarn);
    %--- check suffix
    if length(savename)>4 && ...
         (strcmpi(savename((end-3):end),'_ave') || strcmpi(savename((end-3):end),'-ave') || ...
          strcmpi(savename((end-3):end),'_cov') || strcmpi(savename((end-3):end),'-cov'))
      savename  = savename(1:(end-4));
    end
    cfg.output      = [fullfile(savedir,savename) '_ave' saveext];
    cfg.cov.output  = [fullfile(savedirc,savename) '_cov' saveextc];
  end
  if ~isforce && logical(exist(cfg.cov.output,'file'))      % check overwriting
    curwarn = warning('off', 'backtrace');
    warning('The output file for covariance already exist! Now copying...');
    warning(curwarn);
    copyfile(cfg.cov.output, [cfg.cov.output '.bk']);
  end
end

%-- save output data
if isucorr
    fprintf(1,'Units were corrected...');
end 
fprintf(1,'Writing average...');
fiff_write_evoked(cfg.output, data);
if iscov
    fprintf(1,'covariance...');
    tmpcfg           = cfg.cov;
    tmpcfg.badch     = cfg.badch;
    tmpcfg.unitcorr  = cfg.unitcorr;
    tmpcfg.rawdata   = data.info;
    tmpcfg.overwrite = 'yes';
    tmpcfg.silence   = 'yes';
    ft_mnefif_write_covariance(tmpcfg, varargin{:});
end
fprintf(1,'[done]\n');

path(curpath);
catch ME
    path(curpath);
    rethrow(ME);
end

function I = idx2logical(data, index)
% convert index number into logical vector

I = false(size(data));
I(index) = true;
