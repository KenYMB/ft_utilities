function ft_mnefif_write_covariance(cfg,varargin)
%
% FT_MNEFIF_WRITE_COVARIANCE(cfg,timelock1,timelock2,...) converts fieldtrip
% timelocked average data includeng cov field into mne covariance data, and
% saves fiff file. 
% 
% The configuration should be according to
% 
%   cfg.output     = output file name
%   cfg.badch      = Nx1 cell-array with selection of bad channels (default = ''),
%                    'mag'/'grad'/'vert': select magnetmeter/planar/vertical sensor as bad
%                    '*' is interpreted as a wildcard
%   cfg.merge      = 'yes' or 'no', merge covariance across all input data
%                    correct the degrees of freedom before export (default = 'no')
%   cfg.dofcorr    = 'yes' or 'no', to correct the computational method of
%                    defrees of freedom of covariance, SEE below for further explanation
%                    (default = 'no')
%   cfg.unitcorr   = 'yes' or 'no', to correct magnetometer unit into 'T'
%                    and gradient sensor unit into 'T/m' (default = 'no')
%   *** Caution ***
%   If you load MEG data with ft_preprocessing, the sensor units for
%   planars are considered as T/cm while the intact units for planar
%   sensors are T/m. cfg.unitcorr is only usefull if you corrected this
%   contradiction manually.
% 
%   *** degrees of freedom (dof) ***
%   fieldTrip computes the covariance with the dof of
%       SUM[(nsamples - 1)]_ntrialall
%   MNE computes the covariance with the dof of
%       SUM[nsamples * (ntrialcond - 1)]_ncondition
%   If cfg.doffcorr = 'yes', the covariance is corrected into MNE form.
% 
% Undocumented option
% 
%   cfg.rawdata    = file name of original fiff file (load channel name)
%                    or information structure read by fiff_read_meas_info
%                    'auto': load data.cfg.((previous)).dataset
%                    '':     use data.grad.label or data.label (default)
%   cfg.overwrite  = 'yes' or 'no', (default = 'no')
% 
% ( cfg.diag       = yet not support)
% 
% Exapmle: 
%   cfg                         = [];
%   cfg.preproc.demean          = 'yes';
%   cfg.preproc.baselinewindow  =  [-0.2 0];
%   cfg.covariance              = 'yes';
%   cfg.covariancewindow        = 'all';
%   cfg.removemean              = 'yes';
%   avgMEG = ft_timelockanalysis(cfg,dataMEG);
%   cfg            = [];
%   cfg.output     = 'MEG_cov.fif';
%   cfg.badch      = 'vert';
%   cfg.dofcorr    = 'yes';
%   cfg.unitcorr   = 'yes';                 % reset manual correction
%   ft_mnefif_write_covariance(cfg,avgMEG);
% 
% See also MNE_WRITE_COV, FT_NEUROFIF_WRITE_RAW, FT_NEUROFIF_WRITE_AVERAGE

% using: fieldtrip(, ft_hastoolbox, mne, match_str), cellstrfind, int2ordinal

% 20170829 Yuasa

curpath = path;
ft_defaults;
ft_hastoolbox('mne',1,1);

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
me = 'original:ft_mnefif_write_covariance';
if nargin < 2
    path(curpath);
    error(me,'Incorrect number of arguments');
end

try
%--  prepare cfg
cfg.rawdata    = ft_getopt(cfg,'rawdata','auto');
cfg.output     = ft_getopt(cfg,'output','dummy_cov.fif');
cfg.badch      = ft_getopt(cfg,'badch');
cfg.merge      = ft_getopt(cfg,'merge','no');
cfg.dofcorr    = ft_getopt(cfg,'dofcorr','no');
cfg.unitcorr   = ft_getopt(cfg,'unitcorr','no');
% cfg.diag       = ft_getopt(cfg,'diag',false);
cfg.diag       = false;
cfg.overwrite  = ft_getopt(cfg,'overwrite','no');
cfg.silence    = ft_getopt(cfg,'silence','no');
%--- bad channels
if ~isempty(cfg.badch)                                                  	% interpret cfg.badch
    if  ~iscell(cfg.badch)
        cfg.badch = {cfg.badch};
    else
        cfg.badch = reshape(cfg.badch,1,[]);
    end
    genidx    = match_str(cfg.badch,{'mag','grad','vert'});
    for ich = 1:length(genidx)
        switch cfg.badch{genidx(ich)}
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
%--- unit correction
isucorr  = strcmpi(cfg.unitcorr,'yes');
%--- overwrite flag
isforce  = strcmpi(cfg.overwrite,'yes');
%--- silence flag
issilent = strcmpi(cfg.silence,'yes');
%--- merge & correct covariance
ismerge  = strcmp(cfg.merge,'yes');
isfcorr  = strcmp(cfg.dofcorr,'yes');
ncondinp = nargin - 1;
if ismerge,     ncondout = 1;
else            ncondout = ncondinp;
end

%-- check raw-data
TMPdata = varargin{1};
datasetcand = ft_findcfg(TMPdata.cfg, 'datafile');
if strcmp(cfg.rawdata,'auto')
    if ~isempty(datasetcand);
        cfg.rawdata  = datasetcand;
    else
        warning(me,'Failed to detect raw file.');
        cfg.rawdata  = '';
    end
end

%-- get information of raw-data
datainfo = [];
if ~isempty(cfg.rawdata) && ischar(cfg.rawdata)
  if exist(cfg.rawdata, 'file')
    %--- copy field
    try     datainfo = fiff_read_meas_info(cfg.rawdata);
    catch,  warning(me,'File %s is invalid',cfg.rawdata);
    end
  else
    warning(me,'Cannot open file %s',cfg.rawdata);
  end
end

%-- modify channel info
if isempty(datainfo)
    if isfield(TMPdata,'grad')
        datainfo.ch_names  =  TMPdata.grad.label;
    else
        datainfo.ch_names  =  TMPdata.label;
    end
else
    if isfield(TMPdata,'grad')
        labelmatch = match_str(datainfo.ch_names, TMPdata.grad.label);
    else
        labelmatch = match_str(datainfo.ch_names, TMPdata.label);
    end
    datainfo.ch_names  = datainfo.ch_names(labelmatch);
end
datainfo.ch_names  = reshape(datainfo.ch_names,1,[]);                               % need reshape into row vector
datainfo.nchan     = length(datainfo.ch_names);

%-- time-locked data
cfg.badch        = datainfo.ch_names(cellstrfind(datainfo.ch_names, cfg.badch));    % search & register bad ch
covdat(1:ncondout) = struct('kind',num2cell(1:ncondout),'diag',cfg.diag,'dim',datainfo.nchan,'names',{datainfo.ch_names},'data',0,'projs',struct([]),'bads',{cfg.badch},'nfree',0,'eig',[],'eigvec',[]);
for idats = 1:ncondinp
    %--- check data
    TMPdata = varargin{idats};
    assert(isfield(TMPdata,'cov'),me,'The %s data does not include covariance data', int2ordinal(idats));

    %-- unit correction
    if isucorr && isfield(TMPdata,'grad')
        %--- convert fT into T
        labelmatch     = match_str(TMPdata.grad.label, TMPdata.label);
        corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),{'fT','fT*'});
        TMPdata.cov(corrchan,:)    = TMPdata.cov(corrchan,:) * 1e-15;
        TMPdata.cov(:,corrchan)    = TMPdata.cov(:,corrchan) * 1e-15;
        %--- convert mm/cm/dm into m
        corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),'*/mm');
        TMPdata.cov(corrchan,:)    = TMPdata.cov(corrchan,:) / 1e-3;
        TMPdata.cov(:,corrchan)    = TMPdata.cov(:,corrchan) / 1e-3;
        corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),'*/cm');
        TMPdata.cov(corrchan,:)    = TMPdata.cov(corrchan,:) / 1e-2;
        TMPdata.cov(:,corrchan)    = TMPdata.cov(:,corrchan) / 1e-2;
        corrchan = cellstrfind(TMPdata.grad.chanunit(labelmatch),'*/dm');
        TMPdata.cov(corrchan,:)    = TMPdata.cov(corrchan,:) / 1e-1;
        TMPdata.cov(:,corrchan)    = TMPdata.cov(:,corrchan) / 1e-1;
    end

    %-- nfree computation
    %--   mne:
    %--     for raw data        : (times - 1)
    %--     for epoch data      : ((nave - ncond) * times)
    %--     for keepsamplemean  : (nave * times)
    %--   fieldTrip: 
    %--     cfg.removemean = 'yes': nave * (times - 1)
    %--     cfg.removemean = 'no' : nave * times
    nave  = nanmin(TMPdata.dof(:));
    ntime = ft_findcfg(TMPdata.cfg, 'covariancewindow');
    if isempty(ntime)       % consider as 'all'
        ntime   = length(TMPdata.time);
    else                    % get length of covariance window
        ntime   = nearest(TMPdata.time,ntime);
        ntime   = ntime(2) - ntime(1) + 1;
    end
    isremmean = strcmp(ft_findcfg(TMPdata.cfg, 'removemean'),'yes'); % default = 'yes'
    origfree = nave .* (ntime - isremmean);
    if isfcorr,  TMPfree     = (nave - isremmean) .* ntime; % MNE nfree
    else         TMPfree     = origfree;                    % fieldTrip nfree
    end
    
    %-- 1st step of degrees of freedom correction
    TMPdata.cov         = TMPdata.cov .* origfree;
    
    %-- output cov matrix & nfree
    [labelmatch, labelidx] = match_str(datainfo.ch_names, TMPdata.label);
    if ismerge && idats~=1
      covdat(1).nfree     = covdat(1).nfree + TMPfree;
      if ~covdat(1).diag
         covdat(1).data(labelmatch,labelmatch) = ...
           covdat(1).data(labelmatch,labelmatch) + TMPdata.cov(labelidx,labelidx);
      end
    else
      covdat(idats).nfree = TMPfree;
      if ~covdat(idats).diag
         covdat(idats).data                        = zeros(length(datainfo.ch_names));
         covdat(idats).data(labelmatch,labelmatch) = TMPdata.cov(labelidx,labelidx);
      end
    end
end
%-- 2nd step of degrees of freedom correction
for idats = 1:ncondout
    covdat(idats).data  = covdat(idats).data ./ covdat(idats).nfree;
end

%%%
%-- prepare for output
[savedir,savename,saveext] = fileparts(cfg.output);	% check extention
if ~strcmp(saveext,'.fif') && ~strcmp(saveext,'.fiff')
    cfg.output = [cfg.output '.fif'];
end
if ~isforce && logical(exist(cfg.output,'file'))      % check overwriting
    curwarn = warning('off', 'backtrace');
    warning('The output file for covariance already exist! Now copying...');
    warning(curwarn);
    copyfile(cfg.output, [cfg.output '.bk']);
end

%-- save output data
if ~issilent
 if isucorr
    fprintf(1,'Units were corrected...');
 end 
 fprintf(1,'Writing covariance...');
end
mne_write_cov_all(cfg.output, covdat);
if ~issilent
 fprintf(1,'[done]\n');
end

path(curpath);
catch ME
    path(curpath);
    rethrow(ME);
end

function mne_write_cov_all(fname,covs)
%
%   function mne_write_cov_all(name,covs)
%
%   Write a complete fif file containing a covariance matrix
%
%   fname    filename
%   covs     the covariance matrix to write
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2008/10/10 16:13:57  msh
%   Added mne_ex_read_epochs. Fixed help text of mne_write_cov_file.m
%
%   Revision 1.2  2006/05/03 18:53:06  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.1  2006/04/29 12:44:10  msh
%   Added covariance matrix writing routines.
%

% 20170821 Yuasa: update mne_write_cov_file
%                 fix fclose bug for catch, and enable multiple covs

me='MNE:mne_write_cov_file';

fid = fiff_start_file(fname);

try
    for p = 1:length(covs)
        mne_write_cov(fid,covs(p));
    end
    fiff_end_file(fid);
    return;
catch ME
    fclose(fid);
    delete(fname);
    fprintf('\n');
    rethrow(ME);
end
    
