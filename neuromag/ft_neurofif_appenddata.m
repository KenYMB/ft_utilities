function ft_neurofif_appenddata(cfg)
%
% FT_NEUROFIF_APPENDDATA(cfg)
% 
% cfg.output   = 'filename';    % output file name
% cfg.inputs   = {'filename'};  % load file names for append
% 
% append multiple fif file
% headder information is copied from first fiff file
% 
% % Exapmle: 
% cfg         = [];
% cfg.output  = 'preprocessedMEG.fif';
% cfg.inputs  = {'preprocessedMEG_01.fif';
%                'preprocessedMEG_02.fif'};
% ft_neurofif_appenddata(cfg);

% using: fieldtrip(, ft_hastoolbox, mne, match_str)

% 20161129 Yuasa: create based on mne_ex_read_write_raw.m
% 20161212 Yuasa: minor fix
% 20161219 Yuasa: append bad chs
% 20170112 Yuasa: bug fix
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
    raw = fiff_setup_read_raw(cfg.inputs{1});
catch
    path(curpath);
    error(me,'%s',mne_omit_first_line(lasterr));
end

try
%-- open output files
allch = 1:length(raw.info.ch_names);
[outfid,cals] = fiff_start_writing_raw(cfg.output,raw.info,allch);

%-- append files
for iloop = 1:length(cfg.inputs)
    try
        rawtmp = fiff_setup_read_raw(cfg.inputs{iloop});
    catch
        error(me,'%s',mne_omit_first_line(lasterr));
    end
    %-- channel selection
    [picks, ch_idx] = match_str(raw.info.ch_names, rawtmp.info.ch_names);
    picks = reshape(picks,1,[]);

    %--  Read and write all the data
    for itrl = 1:length(rawtmp.rawdir)
        %-- load raw data
        targetsamp = [rawtmp.rawdir(itrl).first rawtmp.rawdir(itrl).last];  % pickup segment information
        [ refdata, times ] = fiff_read_raw_segment(rawtmp,targetsamp(1),targetsamp(2),allch);
        %-- write new data
        fprintf(1,'Writing...');
        if iloop == 1 && itrl == 1
            fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,targetsamp(1));
        end
        refdata = refdata(ch_idx,:);            % [ch, time]
        fiff_write_raw_buffer(outfid, refdata, cals);
        fprintf(1,'[done]\n');
    end
    if rawtmp.fid > 0
        fclose(rawtmp.fid);
    end
    
    %--  Append bad channels
    badch = match_str(rawtmp.info.bads, raw.info.bads);
    rawtmp.info.bads(badch) = [];
    raw.info.bads = [raw.info.bads rawtmp.info.bads];
end

%--  Apply bad channels
if ~isempty(raw.info.bads)
    badch = match_str(raw.info.ch_names, raw.info.bads);       % sort
    raw.info.bads = raw.info.ch_names(unique(badch));
    fiff_start_block(outfid,FIFF.FIFFB_MNE_BAD_CHANNELS);
    fiff_write_name_list(outfid,FIFF.FIFF_MNE_CH_NAME_LIST,raw.info.bads);
    fiff_end_block(outfid,FIFF.FIFFB_MNE_BAD_CHANNELS);
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