function [grandavg] = ft_connectivitygrandaverage(cfg, varargin)

% FT_CONNECTIVITYGRANDAVERAGE computes the average powerspectrum or time-frequency spectrum
% over multiple subjects
%
% Use as
%   [grandavg] = ft_connectivitygrandaverage(cfg, data1, data2, data3...)
%   [grandavg] = ft_connectivitygrandaverage(cfg, data1_target, data1_control, data2_target,...)     % cfg.diffpair = 'each'
%   [grandavg] = ft_connectivitygrandaverage(cfg, data1_target, data2_target,...,data1_control,...)  % cfg.diffpair = 'group'
%
% The input data is a structure containing the output to FT_CONNECTIVITYANALYSIS
% using a frequency domain metric of connectivity.
% 
% The configuration structure can contain
%   cfg.keepindividual = 'yes' or 'no' (default = 'no')
%   cfg.normalizevar   = 'N' or 'N-1' (default = 'N-1')
%   cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.parameter      = string or cell-array of strings indicating which parameter(s) to average.
%                        default is set to '*spctrm', if it is present in the data.
%   cfg.takeover       = string or cell-array of strings indicating which parameter(s) to take over.
%                        This option only check the first input data.
%   cfg.takediff       = 'yes' or 'no', take difference between paired datas (default = 'no')
%   cfg.diffpair       = 'each' or 'group', how to interpret the datas (default = 'each')
%                        if only one data is input, the subject dimension is considered as different inputs.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input should be
% specified as a cell array.
%
% See also FT_TIMELOCKGRANDAVERAGE, FT_FREQGRANDAVERAGE, FT_CONNECTIVITYANALYSIS

% FIXME averaging coherence is not possible if inputs contain different amounts of data (i.e. chan/freq/time)

% Copyright (C) 2005-2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% Using: cellstrfind

% 20170419 Yuasa: modify FT_FREQGRADAVERAGE
% 20170425 Yuasa: add cfg.takeover
% 20170925 Yuasa: enable takediff for one input data
% 20180713 Yuasa: minor update

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
emplist = false(1,length(varargin));
for i=1:length(varargin)
    if isempty(varargin{i})
        emplist(i) = true;
    else
        varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'no');
    end
end
varargin(emplist) = [];

% set the defaults
cfg.keepindividual = ft_getopt(cfg, 'keepindividual', 'no');
cfg.normalizevar   = ft_getopt(cfg, 'normalizevar', 'N-1');
cfg.channel        = ft_getopt(cfg, 'channel',    'all');
cfg.foilim         = ft_getopt(cfg, 'foilim',     'all');
cfg.toilim         = ft_getopt(cfg, 'toilim',     'all');
cfg.parameter      = ft_getopt(cfg, 'parameter',  '*spctrm');
cfg.takeover       = ft_getopt(cfg, 'takeover',   '');
cfg.takediff       = ft_getopt(cfg, 'takediff',   'no');
cfg.diffpair       = ft_getopt(cfg, 'diffpair',   'each');

fn = fieldnames(varargin{1});
cfg.parameter   = fn(cellstrfind(fn,cfg.parameter));
if isempty(cfg.parameter)
    error('you should specify a valid parameter to average');
end
cfg.takeover   = fn(cellstrfind(fn,cfg.takeover));

if strcmp(cfg.normalizevar,'N')
    cfg.normalizevar = 1;
elseif strcmp(cfg.normalizevar,'N-1')
    cfg.normalizevar = 0;
end

Ndata      = length(varargin);
dimord     = varargin{1}.dimord;
dimtok     = tokenize(dimord,'_');
hasfreq    = ~isempty(strfind(varargin{1}.dimord, 'freq')); % this should always be true
hastime    = ~isempty(strfind(varargin{1}.dimord, 'time'));
hassubj    = ~isempty(strfind(varargin{1}.dimord, 'subj'));
haschancmb = ~isempty(strfind(varargin{1}.dimord, 'chancmb'));
isdiff     = strcmp(cfg.takediff,'yes') & Ndata>1;
singlediff = strcmp(cfg.takediff,'yes') & Ndata==1;


% check whether the input data is suitable
if isdiff&&mod(Ndata,2)
    error('The input datas are not paired');
end

if ischar(cfg.foilim) && strcmp(cfg.foilim, 'all')
    fbeg = -inf;
    fend =  inf;
else
    fbeg = cfg.foilim(1);
    fend = cfg.foilim(2);
end

if ischar(cfg.toilim) && strcmp(cfg.toilim, 'all')
    tbeg = -inf;
    tend =  inf;
else
    tbeg = cfg.toilim(1);
    tend = cfg.toilim(2);
end

% determine which channels, frequencies and latencies are available for all inputs
for i=1:Ndata
    if haschancmb
        varargin{i} = ft_checkdata(varargin{i}, 'cmbrepresentation', 'full');
    end
    cfg.channel = ft_channelselection(cfg.channel, varargin{i}.label);
    if hasfreq
        fbeg = max(fbeg, varargin{i}.freq(1  ));
        fend = min(fend, varargin{i}.freq(end));
    end
    if hastime
        tbeg = max(tbeg, varargin{i}.time(1  ));
        tend = min(tend, varargin{i}.time(end));
    end
end % for i = subject
cfg.foilim = [fbeg fend];
cfg.toilim = [tbeg tend];

% select the data in all inputs
for k=1:numel(cfg.parameter)
    
    % pick the selections
    for i=1:Ndata
        if ~isfield(varargin{i}, cfg.parameter{k})
            error('the field %s is not present in data structure %d', cfg.parameter{k}, i);
        elseif ndims(varargin{i}.(cfg.parameter{k})) ~= length(dimtok)
            error('data.dimord is invalid');
        end
        [dum, chansel] = match_str(cfg.channel, varargin{i}.label);
        varargin{i}.label = varargin{i}.label(chansel);

        if hasfreq
            freqsel = nearest(varargin{i}.freq, fbeg):nearest(varargin{i}.freq, fend);
            varargin{i}.freq = varargin{i}.freq(freqsel);
        end
        if hastime
            timesel = nearest(varargin{i}.time, tbeg):nearest(varargin{i}.time, tend);
            varargin{i}.time = varargin{i}.time(timesel);
        end
        % reshape for average
        if ~hassubj
            varargin{i}.(cfg.parameter{k}) = reshape(varargin{i}.(cfg.parameter{k}),[1 size(varargin{i}.(cfg.parameter{k}))]);
        end
        % select the overlapping samples in the power spectrum
        switch dimord
            case {'chan_freq' 'subj_chan_freq'}
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(:,chansel,freqsel);
            case {'chan_freq_time' 'subj_chan_freq_time'}
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(:,chansel,freqsel,timesel);
            case {'chan_chan_freq' 'chancmb_freq' 'subj_chan_chan_freq' 'subj_chancmb_freq'}
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(:,chansel,chansel,freqsel);
            case {'chan_chan_freq_time' 'chancmb_freq_time' 'subj_chan_chan_freq_time' 'subj_chancmb_freq_time'}
                varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(:,chansel,chansel,freqsel,timesel);
            otherwise
                error('unsupported dimord');
        end
        
    end % for i = subject
end % for k = parameter

% determine the size of the data to be averaged
dim = cell(1,numel(cfg.parameter));
for k=1:numel(cfg.parameter)
    dim{k} = size(varargin{1}.(cfg.parameter{k}));
    dim{k}(1) = [];     % remove number of subj
end

% compute Nsubj
if isdiff
    eachsubj = zeros(1,Ndata/2);
    for i=2:2:Ndata
      eachsubj(i/2) = size(varargin{i}.(cfg.parameter{1}),1);
      if eachsubj(i/2) ~= size(varargin{i-1}.(cfg.parameter{1}),1)
          error('The input datas are not paired');
      end
    end
    Nsubj = sum(eachsubj);
elseif singlediff
    eachsubj = size(varargin{1}.(cfg.parameter{1}),1);
    if mod(eachsubj,2)
        error('The input datas are not paired');
    else
        eachsubj = eachsubj/2;
    end
    Nsubj = eachsubj;
else
    eachsubj = zeros(1,Ndata);
    for i=1:Ndata
      eachsubj(i) = size(varargin{i}.(cfg.parameter{1}),1);
    end
    Nsubj = sum(eachsubj);
end

% give some feedback on the screen
if strcmp(cfg.keepindividual, 'no')
    for k=1:numel(cfg.parameter)
        fprintf('computing average %s over %d subjects\n', cfg.parameter{k}, Nsubj);
    end
else
    for k=1:numel(cfg.parameter)
        fprintf('not computing average, but keeping individual %s for %d subjects\n', cfg.parameter{k}, Nsubj);
    end
end

% allocate memory to hold the data and collect it
for k=1:numel(cfg.parameter)
    % marge subjects
    tmp = zeros([Nsubj dim{k}]);
    if isdiff
      switch cfg.diffpair
        case 'each'
          for s=2:2:Ndata
            tmp((sum(eachsubj(1:(s/2-1)))+1):sum(eachsubj(1:(s/2))),:,:,:,:) =...
                varargin{s-1}.(cfg.parameter{k}) - varargin{s}.(cfg.parameter{k});
          end
        case 'group'
          for s=1:(Ndata/2)
            tmp((sum(eachsubj(1:(s-1)))+1):sum(eachsubj(1:s)),:,:,:,:) =...
                varargin{s}.(cfg.parameter{k}) - varargin{s+Ndata/2}.(cfg.parameter{k});
          end
        otherwise
          error('invalid parameter of cfg.diffpair');
      end
    elseif singlediff
      switch cfg.diffpair
        case 'each'
            tmp = varargin{1}.(cfg.parameter{k})(1:Nsubj,:,:,:,:) ...
                - varargin{1}.(cfg.parameter{k})(2:Nsubj,:,:,:,:);
        case 'group'
            tmp = varargin{1}.(cfg.parameter{k})(1:Nsubj,:,:,:,:) ...
                - varargin{1}.(cfg.parameter{k})((Nsubj+1):end,:,:,:,:);
        otherwise
          error('invalid parameter of cfg.diffpair');
      end
    else
      for s=1:Ndata
        tmp((sum(eachsubj(1:(s-1)))+1):sum(eachsubj(1:s)),:,:,:,:) =...
            varargin{s}.(cfg.parameter{k});
      end
    end
    % take average
    if strcmp(cfg.keepindividual, 'no')
        grandavg.(cfg.parameter{k})          = reshape(nanmean(tmp,1),dim{k});
        grandavg.([cfg.parameter{k} '_var']) = reshape(nanvar(tmp,cfg.normalizevar,1),dim{k});
    else
        grandavg.(cfg.parameter{k}) = tmp;
    end
end

% collect the output data
grandavg.label = varargin{1}.label;
grandavg.freq  = varargin{1}.freq;
if isfield(varargin{1}, 'time')
    % remember the time axis
    grandavg.time = varargin{1}.time;
end
if isfield(varargin{1}, 'grad')
    warning('discarding gradiometer information because it cannot be averaged');
end
if isfield(varargin{1}, 'elec')
    warning('discarding electrode information because it cannot be averaged');
end
if strcmp(cfg.keepindividual, 'yes')
    grandavg.dimord = ['subj_',varargin{1}.dimord];
elseif strcmp(cfg.keepindividual, 'no')
    cfg.subjects  = Nsubj;
    if hassubj
        grandavg.dimord = varargin{1}.dimord(6:end);
    else
        grandavg.dimord = varargin{1}.dimord;
    end
end
if haschancmb
    grandavg = ft_checkdata(grandavg, 'cmbrepresentation', 'sparse');
end
for k=1:numel(cfg.takeover)
    grandavg.(cfg.takeover{k})  = varargin{1}.(cfg.takeover{k});
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance grandavg
ft_postamble history    grandavg
ft_postamble savevar    grandavg
