function ft_neurofif_write_raw(cfg,data)
%
% FT_NEUROFIF_WRITE_RAW(cfg,data)
%   Output fieldtrip data as fif file
%   NaN or 0-filled trials are removed as bad trial
%   NaN or 0-filled channels are filled 0 and marked as bad channel
%
% data      % fieldtrip data
%           % NaN trials, NaN channels are skipped
% 
% cfg.output   = 'filename';    % output file name
% cfg.input    = 'auto';        % loaded file name (header information of output refers this file)
%                               % 'auto': load data.cfg.[previous].dataset
% cfg.trig     = 'yes';         % add triger data : please prepare data.trig = {trials}chs*times
% cfg.gtrig    = 'yes';         % add general triger data : please prepare data.gtrig = {trials}chs*times
% cfg.badch    = {'ch_name'};   % marked as bad channels (don't take over raw bad ch information)
%                               % 'mag'/'grad'/'vert': select magnetmeter/planar/vertical sensor as bad
%                               % '*' is interpreted as a wildcard
% 
% -- unusual option --
% cfg.sort     = 'no';          % if 'yes', sort trial numbers based on sampleinfo
% cfg.appended = 'auto';        % if 'yes', close gaps of sampleinfo (unusual case)
%                               % this option is usefull for appended data which contains bad trials
%                               % if 'no', use sampleinfo to reconstruct fiff file
%                               % if 'auto', judge if data is appended based on data.cfg.previous
% cfg.copysegment = 'no';       % if 'yes', skip the warning if the raw data length is less than MATLAB data
%                               % (automatically answer 'y' & fill dummy data)
% cfg.unitcorr    = 'no';       % if 'yes', correct magnetometer unit into 'T'
%                               % and gradient sensor unit into 'T/m' based on data.grad
%                               %-------- Caution -------
%                               % If you load MEG data with ft_preprocessing, the sensor units for
%                               % planars are considered as T/cm while the intact units for planar
%                               % sensors are T/m. cfg.unitcorr is only usefull if you corrected this
%                               % contradiction manually.
% 
% -- Undocumented option --
% cfg.overwrite  = 'no'
% 
% 
% Exapmle: 
%   dataMEG.trig  = trigger_data;
%   dataMEG.gtrig = general_trigger_data;
%   cfg           = [];
%   cfg.output    = 'preprocessedMEG.fif';
%   cfg.badch     = 'vert';
%   ft_neurofif_write_raw(cfg,dataMEG);
% 
% See also FT_NEUROFIF_READ_RAW, FT_NEUROFIF_WRITE_AVERAGE

% using: fieldtrip(, ft_hastoolbox, mne, match_str), cellstrfind

% 20160824 Yuasa: create based on mne_ex_read_write_raw.m
% 20160830 Yuasa: sort trials based on sampleinfo
% 20161024 Yuasa: add general trigger
% 20161129 Yuasa: disable dimord warning in ft_preprocessing
% 20161212 Yuasa: bug fix (trig & gtrig), change bad chs processing
% 20161219 Yuasa: enable wildcard for bad chs
% 20170112 Yuasa: correct data.sampleinfo for appended data (add cfg.appended)
% 20170113 Yuasa: error for the case: matlab data size > load data size
% 20170126 Yuasa: add option: ignore raw data length & reuse 1st segment
% 20170217 Yuasa: help update
% 20170623 Yuasa: minor update
% 20170821 Yuasa: implement cfg.unitcorr option
% 20180501 Yuasa: bug fix for auto cfg.input detection
% 20190508 Yuasa: enable '-' option in cfg.badch


curpath = path;
ft_defaults;
ft_hastoolbox('mne',1,1);

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
me = 'original:ft_neurofif_write_raw';
if nargin ~= 2
    path(curpath);
    error(me,'Incorrect number of arguments');
end

try
%--  prepare cfg
issort           = isfield(cfg,'sort') && strcmpi(cfg.sort,'yes');    % default = 'no'
istrigger        = ~isfield(cfg,'trig') || ~strcmpi(cfg.trig,'no');   % default = 'yes'
isgentrigger     = ~isfield(cfg,'gtrig') || ~strcmpi(cfg.gtrig,'no'); % default = 'yes'
isbadch          = isfield(cfg,'badch') && ~isempty(cfg.badch) && ...
                    (iscellstr(cfg.badch) || ischar(cfg.badch));      % default = []
iscopysegment    = isfield(cfg,'copysegment') && ...
                    strcmpi(cfg.copysegment,'yes');                   % default = 'no'(hidden parameter)
isucorr          = isfield(cfg,'unitcorr') && strcmpi(cfg.unitcorr,'yes');   % default = 'no';
isforce          = isfield(cfg,'overwrite') && strcmpi(cfg.overwrite,'yes'); % default = 'no';
if istrigger,              trigdata = data.trig;    end
if isgentrigger,           gtrigdata = data.gtrig;  end
if ~isfield(cfg,'input'),  cfg.input    = 'auto';   end               % default = 'auto'
if ~isfield(cfg,'appended'),  cfg.appended    = 'auto';   end         % default = 'auto'
isappended       = false;
datasetcand      = [];
if isfield(data,'cfg')          % SEE ALSO ft_findcfg(data.cfg, 'datafile')
    precfg = data.cfg;
    while isfield(precfg,'previous')
        if iscell(precfg.previous) && ~isempty(precfg.previous)
            precfg = precfg.previous{1};
            isappended = length(precfg.previous) ~= 1;
        elseif ~iscell(precfg.previous) && ~isempty(precfg.previous)
            precfg = precfg.previous;
        else
            precfg = rmfield(precfg, 'previous');
        end
        if isfield(precfg,'dataset')
            datasetcand  = precfg.dataset;
        end
    end
end
if strcmp(cfg.input,'auto')
    if isfield(precfg,'dataset')
        cfg.input  = precfg.dataset;
    elseif ~isempty(datasetcand);
        cfg.input  = datasetcand;
    else
        error(me,'Failed to detect raw file. Please specify cfg.input.');
    end
end
if strcmp(cfg.appended,'yes')
    isappended      = true;
elseif strcmp(cfg.appended,'no')
    isappended      = false;
end %-- skip when cfg.appended = 'auto' (check in above part)

%-- check overwrite
if ~isforce && logical(exist(cfg.output,'file'))
    curwarn = warning('off', 'backtrace');
    warning('The output file already exist! Now copying...');
    warning(curwarn);
    copyfile(cfg.output, [cfg.output '.bk']);
end

%-- reject bad trials
    orig_state = warning;
    warning('off','all');
tmpcfg                 = [];
tmpcfg.trials          = squeeze(nanmax(nanmax(abs(cat(3,data.trial{:})),[],1),[],2));  % NaN / 0 trial detection
tmpcfg.trials          = ~logical(isnan(tmpcfg.trials) + (tmpcfg.trials==0));           % [NaN 0] trial rejection
data                   = ft_preprocessing(tmpcfg,data);
if isappended                                            % close gaps of sampleinfo for appended data
    data               = rmfield(data,'sampleinfo');
    data               = ft_preprocessing([],data);
end 
    warning(orig_state);
if istrigger,          trigdata   = trigdata(tmpcfg.trials);    end
if isgentrigger,       gtrigdata  = gtrigdata(tmpcfg.trials);   end

%--  Setup for reading the raw data
try
    raw = fiff_setup_read_raw(cfg.input);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
%-- error judgement
lastsample  = max(data.sampleinfo(:));
if lastsample > double(raw.last_samp - raw.first_samp + 1)
    warning(me, 'Loaded data size is smaller than input data size.');
    if ~iscopysegment
        iscopysegment = input('Do you want to copy 1st segment information of loaded data for all trials? (yes/no) :', 's');
        fprintf(1,'You can use cfg.copysegment\n');
        iscopysegment = logical(max(strcmpi(iscopysegment,{'yes';'y';'1'})));
    end
    assert(iscopysegment, me, 'Prepare raw data which contains all epochs you want to write.');
else
    iscopysegment = false;
end

%-- channel selection
[picks, ch_idx] = match_str(raw.info.ch_names, data.label);
picks = reshape(picks,1,[]);
allch = 1:length(raw.info.ch_names);
%-- trigger channel
if istrigger
    picksadd = reshape(cellstrfind(raw.info.ch_names, 'STI0*'),1,[]);
    if length(picksadd) < size(trigdata{1},1)
        error(me,'too match triggers');
    end
    picksadd = picksadd(1:size(trigdata{1},1));
    picks = [picks picksadd];
end
%-- general trigger channel
if isgentrigger
    picksadd = reshape(cellstrfind(raw.info.ch_names, 'STI1*'),1,[]);
    if length(picksadd) < size(gtrigdata{1},1)
        error(me,'too match triggers');
    end
    picksadd = picksadd(1:size(gtrigdata{1},1));
    picks = [picks picksadd];
end
%-- bad channels
nanchs         = ~squeeze(nanmax(nanmax(abs(cat(3,data.trial{:})),[],3),[],2));   % NaN / 0 ch detection
nanchs         = logical(isnan(nanchs) + (nanchs));                               % [NaN 0] ch rejection
for itrl = 1:length(data.trial)
    data.trial{itrl}(nanchs,:) = 0;                                               % 0-fill for NaN / 0 ch
end
nanchs         = cellstrfind(raw.info.ch_names,data.label(nanchs));              % sort based on header info
    nanchs     = idx2logical(raw.info.ch_names,nanchs);
if isbadch                                                                       % add indicated bad ch
    if ischar(cfg.badch)
        cfg.badch = {cfg.badch};
    else
        cfg.badch = reshape(cfg.badch,1,[]);
    end
    genidx  = cellstrfind(cfg.badch,{'all','mag','grad','vert'});
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
    %-- exclude '-*' channels
    negidx = ~cellfun(@isempty,regexp(cfg.badch,'^-','once'));
    negchs = cellstrfind(raw.info.ch_names, regexprep(cfg.badch(negidx),'^-',''));
        negchs     = idx2logical(raw.info.ch_names,negchs);
    badchs = cellstrfind(raw.info.ch_names, cfg.badch(~negidx));
        badchs     = idx2logical(raw.info.ch_names,badchs);
    badchs  = badchs & ~negchs;                                                    % exclude '-*' chs from bad-channels
    cfg.badch = raw.info.ch_names(badchs);
    
end
raw.info.bads    = raw.info.ch_names(badchs | nanchs);                             % register bad ch

%-- unit correction
if isucorr
  if isfield(data,'grad')
     convval = ones(size(data.trial{1}));
     %--- convert fT into T
     labelmatch     = match_str(data.grad.label, data.label);
     corrchan = cellstrfind(data.grad.chanunit(labelmatch),{'fT','fT*'});
     convval(corrchan,:)     = convval(corrchan,:) * 1e-15;
     %--- convert mm/cm/dm into m
     corrchan = cellstrfind(data.grad.chanunit(labelmatch),'*/mm');
     convval(corrchan,:)     = convval(corrchan,:) / 1e-3;
     corrchan = cellstrfind(data.grad.chanunit(labelmatch),'*/cm');
     convval(corrchan,:)     = convval(corrchan,:) / 1e-2;
     corrchan = cellstrfind(data.grad.chanunit(labelmatch),'*/dm');
     convval(corrchan,:)     = convval(corrchan,:) / 1e-1;
  else
     isucorr = false;
     warning(me, 'Cannot specify current units, ''grad'' field does not exist');
  end
end

%-- sort info
if issort
    trialseqnum = [data.sampleinfo [1:size(data.sampleinfo,1)]'];
    trialseqnum = sortrows(trialseqnum,1);
    trialseqnum = trialseqnum(:,end);
else
    trialseqnum = [1:size(data.sampleinfo,1)]';
end

%-- open output files
[outfid,cals] = fiff_start_writing_raw(cfg.output,raw.info,allch);

%--  Read and write all the data
for itrl = 1:length(data.trial)
    %-- load raw data
    if iscopysegment
      targetsamp = int32(data.sampleinfo(trialseqnum(1),:) + double(raw.first_samp) - 1);  % 1st segment
    else
      targetsamp = int32(data.sampleinfo(trialseqnum(itrl),:) + double(raw.first_samp) - 1);  % correct sample info between fieldtrip and neuromag
    end
    [ refdata, times ] = fiff_read_raw_segment(raw,targetsamp(1),targetsamp(2),allch);
    %--- apply unit correction
    if isucorr
      fprintf(1,'Correcting units...');
      data.trial{trialseqnum(itrl)}   =  data.trial{trialseqnum(itrl)} .* convval;
    end
    %-- write new data
    fprintf(1,'Writing...');
    if trialseqnum(itrl) == 1
        fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,targetsamp(1));
    elseif (trialseqnum(itrl) == length(data.trial)) && ~isappended && ~iscopysegment
        fiff_write_int(outfid,FIFF.FIFF_LAST_SAMPLE,targetsamp(2));
    end
    rplcdata = data.trial{trialseqnum(itrl)}(ch_idx,:);
    if istrigger, 
        rplcdata = cat(1,rplcdata,trigdata{trialseqnum(itrl)});     end
    if isgentrigger, 
        rplcdata = cat(1,rplcdata,gtrigdata{trialseqnum(itrl)});    end
    refdata(picks,:) = rplcdata;
    fiff_write_raw_buffer(outfid, refdata, cals);
    fprintf(1,'[done]\n');
end

fiff_finish_writing_raw(outfid);
if raw.fid > 0
    fclose(raw.fid);
end

path(curpath);
catch ME
    path(curpath);
    rethrow(ME);
end

function I = idx2logical(data, index)
% convert index number into logical vector

I = false(size(data));
I(index) = true;
