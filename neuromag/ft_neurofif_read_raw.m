function data = ft_neurofif_read_raw(cfg)
%
% data = FT_NEUROFIF_READ_RAW(cfg)
% 
%  data      = fieldtrip data
%  
%  The congfigurations should be specified with
%   cfg.dataset      = string with the filename
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.trigchan     = string with the trigger channel name (default = 'STI*')
%                      end of trigchan is available if multiple channels are input
%   cfg.keeptrig     = 'no' or 'yes'  return trial data of trigger channel or not (default = 'no')
%   cfg.offset       = number or Nx1 array with offset time for each trial [s]
% 
% 
% Read fif file which was created by FT_NEUROFIF_WRITE_RAW
% in the form of fieldtrip
% 
% Exapmle: 
%   cfg           = [];
%   cfg.output    = 'preprocessedMEG.fif';
%   cfg.channel   = 'MEG';
%   cfg.offset    = 0.2;                      % start at -0.2s
%   dataMEG = ft_neurofif_write_raw(cfg);
% 
% See also FT_NEUROFIF_WRITE_RAW

% using: fieldtrip(, ft_hastoolbox, mne, match_str)

% 20170623 Yuasa: create based on mne_ex_read_write_raw.m
% 20170810 Yuasa: minor update


%-- these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

%-- do the general setup of the function
curpath = path;
ft_defaults
ft_hastoolbox('mne',1,1);
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

%-- the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

%-- neuromag setting
global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
me = 'original:ft_neurofif_write_raw';

try
%-- set the defaults
cfg.channel       = ft_getopt(cfg, 'channel', 'all');      
cfg.trigchan      = ft_getopt(cfg, 'trigchan', 'STI*');
cfg.offset        = ft_getopt(cfg, 'offset', 0);
cfg.keeptrig      = ft_getopt(cfg, 'keeptrig', 'no');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');        % is passed to low-level function, empty implies autodetection
cfg.dataformat    = ft_getopt(cfg, 'dataformat');          % is passed to low-level function, empty implies autodetection
cfg.coordsys      = ft_getopt(cfg, 'coordsys', 'head');    % is passed to low-level function
cfg.coilaccuracy  = ft_getopt(cfg, 'coilaccuracy');        % is passed to low-level function

%-- check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required',   {'headerfile', 'datafile'});
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

%-- read the header
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'coordsys', cfg.coordsys, 'coilaccuracy', cfg.coilaccuracy);
rawinfo   = hdr.orig.raw;   % SEE ALSO, fiff_setup_read_raw(cfg.dataset)

%-- prepare trial info
ntrl  = length(rawinfo.rawdir);
Fs       = hdr.Fs;
if numel(cfg.offset)==1
    cfg.offset = repmat(cfg.offset,1,ntrl);
end
assert(numel(cfg.offset)==ntrl, 'The structure of ''cfg.offset'' is invalid.');

%-- prepare trial
trl      = zeros(ntrl, 3);
for itrl=1:ntrl
    trl(itrl,1) = rawinfo.rawdir(itrl).first - rawinfo.first_samp + 1;
    trl(itrl,2) = rawinfo.rawdir(itrl).last  - rawinfo.first_samp + 1;
    trl(itrl,3) = round(-cfg.offset(itrl).*Fs);
end

%-- prepare channel
cfg.channel  = ft_channelselection(cfg.channel, hdr.label);
cfg.trigchan = ft_channelselection(cfg.trigchan, hdr.label);
[chnidx, rawchn] = match_str(cfg.channel, hdr.label);
[trgidx, rawtrg] = match_str(cfg.trigchan, hdr.label);
nchn = length(chnidx);  ntrg = length(trgidx);

%-- main
cutdata     = cell(1, ntrl);
timedata    = cell(1, ntrl);
trigdata    = cell(1, ntrl);
triginfo    = ones(1,length(rawinfo.rawdir));               % fill 1, because conditions are already segrigated
for itrl = 1:length(rawinfo.rawdir)
    %-- load raw data
    targetsamp = [rawinfo.rawdir(itrl).first rawinfo.rawdir(itrl).last];  % pickup segment information
    [refdata, times] = fiff_read_raw_segment(rawinfo,targetsamp(1),targetsamp(2),[rawchn; rawtrg]);
    cutdata{itrl}    = refdata(1:nchn,:);
    trigdata{itrl}   = refdata((1:ntrg)+nchn,:);
    triginfo(itrl)   = max(trigdata{itrl}(end,:));
    timedata{itrl}   = times - times(1) + cfg.offset(itrl);
end

%-- output
if isfield(hdr, 'orig')
    hdr = rmfield(hdr, 'orig');
end
cfg.trl     = trl;

dataout                    = [];
dataout.hdr                = hdr;                  % header details of the datafile
dataout.label              = hdr.label(rawchn);    % labels of channels that have been read, can be different from labels in file due to montage
dataout.time               = timedata;             % vector with the timeaxis for each individual trial
dataout.trial              = cutdata;
if strcmpi(cfg.keeptrig,'yes')
   dataout.trigtrial       =  trigdata;
end
dataout.fsample            = Fs;
dataout.sampleinfo         = cfg.trl(:,1:2);
if size(cfg.trl,2) > 3
  dataout.trialinfo        = cfg.trl(:,4:end);
end
dataout.triginfo           = triginfo;
if isfield(hdr, 'grad')
  dataout.grad             = hdr.grad;             % gradiometer system in head coordinates
end
if isfield(hdr, 'elec')
  dataout.elec             = hdr.elec;             % EEG information in header (f.e. headerformat = 'neuromag_fif')
end
if isfield(hdr, 'opto')
  dataout.opto             = hdr.opto;             % NIRS  information in header (f.e. headerformat = 'artinis')
end

%-- do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

%-- rename the output variable to accomodate the savevar postamble
data = dataout;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

path(curpath);
catch ME
    %-- do the general cleanup and bookkeeping at the end of the function
    ft_postamble debug
    ft_postamble trackconfig
    
    path(curpath);
    rethrow(ME);
end