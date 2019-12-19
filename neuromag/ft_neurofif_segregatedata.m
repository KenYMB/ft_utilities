
function ft_neurofif_segregatedata(cfg)
%
% FT_NEUROFIF_SEGREGATEDATA(cfg)
% 
% cfg.outputs   = {'filename'};     % output file name
%                                   % 'auto': add suffix of '_01','_02',...
% cfg.input     = 'filename';       % load file names for append
% cfg.sessions  = [startpoint endpoint];
%                   % [n x 2] matrix: indicate by timepoints
% 
% segregate fif file into multiple fif file
% output files are segrigated at the border of epochs
% 
% % Exapmle: 
% cfg         = [];
% cfg.outputs = {'preprocessedMEG_01.fif';
%                'preprocessedMEG_02.fif'};
% cfg.input   = 'preprocessedMEG_0102.fif';
% cfg.sessions  = [1 50000; 50001 100000];
% ft_neurofif_segregatedata(cfg);

% using: fieldtrip(, ft_hastoolbox, mne), nearlyeq

% 20161130 Yuasa: create based on mne_ex_read_write_raw.m
% 20161212 Yuasa: minor fix
% 20170113 Yuasa: fix for the case: indicated endpoint > data size
% 20170810 Yuasa: minor update

curpath = path;
ft_defaults;
ft_hastoolbox('mne',1,0);

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
me = 'original:ft_neurofif_appenddata';
if nargin ~= 1
    path(curpath);
    error(me,'Incorrect number of arguments');
end

%--  Setup for reading the raw data
try
    raw = fiff_setup_read_raw(cfg.input);
catch
    path(curpath);
    error(me,'%s',mne_omit_first_line(lasterr));
end

%-- fix timepoints info
cfg.sessions                    = int32(cfg.sessions) + int32(raw.first_samp);
if max(cfg.sessions(:) < int32(raw.first_samp)) || max(cfg.sessions(:) > int32(raw.last_samp))
    cfg.sessions(cfg.sessions < int32(raw.first_samp)) = int32(raw.first_samp);
    cfg.sessions(cfg.sessions > int32(raw.last_samp))  = int32(raw.last_samp);
    warning(me,'cfg.sessions is corrected');
end

try
%-- Prepare output files
if ~isfield(cfg,'outputs') || (ischar(cfg.outputs) && strcmp(cfg.outputs,'auto'))
    [basedir, basename, baseext] = fileparts(cfg.input);
elseif ischar(cfg.outputs)        % ~iscellstr(cfg.outputs)
    [basedir, basename, baseext] = fileparts(cfg.outputs);
end
if ~isfield(cfg,'outputs') || ischar(cfg.outputs)
    if ~isempty(basedir), basename = [basedir '/' basename];    end
    cfg.outputs = cell(size(cfg.sessions,1),1);
    for iloop = 1:size(cfg.sessions,1)
        cfg.outputs{iloop} = sprintf('%s_%02d.fif',basename,iloop);
    end
end

assert(size(cfg.sessions,1)==length(cfg.outputs),me,'Prepare output files for each session!');

%-- channel selection
picks = raw.info.ch_names;
allch = 1:length(raw.info.ch_names);

%-- append files
for iloop = 1:size(cfg.sessions,1)
    fprintf(1,'Start to write ''%s''\n',cfg.outputs{iloop});
    %-- open output files
    [outfid,cals] = fiff_start_writing_raw(cfg.outputs{iloop},raw.info,allch);
    
    %-- trial selection
    trials = [nearlyeq([raw.rawdir(:).first], cfg.sessions(iloop,1),'small',1), ...
              nearlyeq([raw.rawdir(:).last], cfg.sessions(iloop,2),'large',1)];
    %--  Read and write all the data
    for itrl = trials(1):trials(2)
        %-- load raw data
        targetsamp = [raw.rawdir(itrl).first raw.rawdir(itrl).last];  % pickup segment information
        [ refdata, times ] = fiff_read_raw_segment(raw,targetsamp(1),targetsamp(2),allch);
        %-- write new data
        fprintf(1,'Writing...');
        if itrl == trials(1)
            fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,targetsamp(1));
        end
        fiff_write_raw_buffer(outfid, refdata, cals);
        fprintf(1,'[done]\n');
    end
    
    fiff_finish_writing_raw(outfid);
end

if raw.fid > 0
    fclose(raw.fid);
end

path(curpath);
catch ME
    path(curpath);
    rethrow(ME);
end