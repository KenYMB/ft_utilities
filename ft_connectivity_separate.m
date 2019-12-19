function      [data] = ft_connectivity_separate(data,separatetype,varargin)
% FT_CONNECTIVITY_SEPARATE pickup feedforward/feedback channels 
% 
% Usage:
%   [data] = ft_connectivity_separate(data,'feedforward')
%   [data] = ft_connectivity_separate(data,'feedback','cmbrepresentation', 'sparse')
% 
% Options: 
%   1st option         = feedforward, feedback, both, all
%   cmbrepresentation  = sparse, full (reshape matrix of covariance and cross-spectral density)
%                        default: estimate and keep current data representation 
%   channel            = all, cell-array (see also ft_channelselection)
%                        default: all
%   keepchannel        = yes, no (if keep unselected channels)
%                        default: yes

% 20170424 Yuasa
% 20180517 Yuasa: add error for typo; add 'keepchannel'

%-- these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

%-- do the general setup of the function
ft_defaults

%-- get the optional input arguments
narginchk(2,inf);
channel              = ft_getopt(varargin, 'channel','all');
keepchannel          = ft_getopt(varargin, 'keepchannel','yes');
cmbrepresentation    = ft_getopt(varargin, 'cmbrepresentation','');
assert(any(strcmp(separatetype,{'feedforward','feedback','both','all'})),'specified 1st option is invalid')
assert(any(strcmp(cmbrepresentation,{'sparse','full',''})),'cmbrepresentation must be ''sparse'' or ''full''')

%-- check dimord
dimord      = data.dimord;
dimtok      = tokenize(dimord, '_');
dim_chan    = find(strcmp(dimtok,'chan'));
dim_chancmb = find(strcmp(dimtok,'chancmb'));
ndimdat     = length(dimtok);
if any(strcmp('chan',  dimtok)),    nchn=length(data.label);       else nchn = 0;    end
if any(strcmp('chancmb',  dimtok)), nchncmb=length(data.labelcmb); else nchncmb = 0; end
if any(strcmp('freq',  dimtok)),    nfrq=length(data.freq);        else nfrq = 1;    end
if any(strcmp('time',  dimtok)),    ntim=length(data.time);        else ntim = 1;    end

if length(dim_chan)~=2 && isempty(dim_chancmb)
    error('the data is not a connectivity data');
end

%-- specify data.parameter
fn = fieldnames(data);
for ilp =length(fn):-1:1
  okflag = true;
  if ndims(data.(fn{ilp})) ~= ndimdat, 	okflag = false;
  else
    for idims =1:ndimdat
      switch dimtok{idims}
         case 'chan'
           if size(data.(fn{ilp}),idims) ~= nchn,  okflag = false;    end
         case 'chancmb'
           if size(data.(fn{ilp}),idims) ~= nchncmb,  okflag = false; end
         case 'freq'
           if size(data.(fn{ilp}),idims) ~= nfrq,  okflag = false;    end
         case 'time'
           if size(data.(fn{ilp}),idims) ~= ntim,  okflag = false;    end
      end
    end
  end
  if ~okflag,   fn(ilp) = [];   end
end
if isempty(fn)
    error('data.dimord is invalid');
end

%-- channel
if iscellstr(channel)
  nominus = false;
  for ilp=1:length(channel)
     if ~strcmp(channel{ilp}(1),'-'),   nominus = true; end
  end
  if ~nominus,  channel = [{'all'}; channel(:)];    end
end

%-- prepare
if ~isempty(dim_chancmb)
    if isempty(cmbrepresentation),  cmbrepresentation = 'sparse';   end
    %-- rearrange labels
    if isfield(data,'label')
        deslabel  = data.label;
    else
        deslabel  = unique(data.labelcmb,'stable');
        findfirst = match_str(data.labelcmb(:,1),deslabel(1));
        if ~strcmp(data.labelcmb{findfirst(1),:})       % if {label1, label1} is already removed
            deslabel = [deslabel(end); deslabel(1:(end-1))];
        end
        data.label   = deslabel;
    end
    nchn = length(data.label);
    [~,labellist]  = match_str(data.labelcmb,deslabel);
    labellist = [reshape(labellist,[],2), (1:size(labellist)/2)'];
    for ichn = 1:nchn
       for jchn = 1:nchn
          if ~any(and(labellist(:,1)==ichn, labellist(:,2)==jchn))
              labellist = [labellist; ichn, jchn, nchncmb+1];
          end
       end
    end
    labellist = sortrows(labellist,[2 1]);
    %-- set nan
    strforpick = '(';
    for idm =1:ndimdat;
        if idm~=1,  strforpick = [strforpick, ', '];    end
        if idm == dim_chancmb
            strforpick = [strforpick, 'nchncmb+1'];
        else
            strforpick = [strforpick, ':'];
        end
    end
    strforpick = [strforpick, ')'];
    for ilp =1:length(fn)
       if islogical(data.(fn{ilp}))
         eval(sprintf('data.(fn{ilp})%s = false;',strforpick));
       else
         eval(sprintf('data.(fn{ilp})%s = nan;',strforpick));
       end
    end
    %-- set full label space
    strforpick = '(';
    for idm =1:ndimdat;
        if idm~=1,  strforpick = [strforpick, ', '];    end
        if idm == dim_chancmb
            strforpick = [strforpick, 'labellist(:,3)'];
        else
            strforpick = [strforpick, ':'];
        end
    end
    strforpick = [strforpick, ')'];
    for ilp =1:length(fn)
       eval(sprintf('data.(fn{ilp}) = data.(fn{ilp})%s;',strforpick));
    end
    
     sizparam = size(data.(fn{1}));
     nchncmb = nchn^2;
else
   %-- reshape data
   if isempty(cmbrepresentation),  cmbrepresentation = 'full';   end
   nchncmb = nchn^2;     dim_chancmb = dim_chan(1);
   ndimdat = ndimdat - 1;
   for ilp =1:length(fn)
     sizparam = size(data.(fn{ilp}));
     sizparam(dim_chan(2)) = [];  sizparam(dim_chan(1)) = nchncmb;
     data.(fn{ilp}) = reshape(data.(fn{ilp}), sizparam);
   end
   dimtok(dim_chan(2)) = [];  dimtok{dim_chan(1)} = 'chancmb';
end
data.labelcmb = [reshape(repmat(reshape(data.label,[],1),1,nchn),[],1),...
                 reshape(repmat(reshape(data.label,1,[]),nchn,1),[],1)];

%-- get index
switch separatetype
    case 'feedforward'
        pickidx = find(~triu(ones(nchn,nchn),1));       % tilde for nan
    case 'feedback'
        pickidx = find(~tril(ones(nchn,nchn),-1));      % tilde for nan
    case 'all'
        pickidx = find(diag(ones(1,nchn)));             % not tilde for nan
    case 'both'
        data.dimord     = sprintf('%s_',dimtok{:});  data.dimord(end)=[];
        optflg = {'channel',channel,'cmbrepresentation',cmbrepresentation};
        tempdata = ft_connectivity_separate(data,'feedforward',optflg{:});      % feedforward
        switch cmbrepresentation
          case 'sparse'
            fullsizdata = [sizparam(1:(dim_chancmb-1)), nchn, nchn, sizparam((dim_chancmb+1):end)];
            fulldimperm = 1:(ndimdat+1);    fulldimperm([0 1]+dim_chancmb) = fulldimperm([1 0]+dim_chancmb);
            for ilp=1:length(fn)    % replace feedforward by feedback
                data.(fn{ilp}) = reshape(permute(reshape(data.(fn{ilp}),fullsizdata),fulldimperm),sizparam);
            end
            data     = ft_connectivity_separate(data,'feedforward',optflg{:});  % feedback
          case 'full'
            data     = ft_connectivity_separate(data,'feedback',optflg{:});     % feedback
        end
        for ilp =1:length(fn)
            data.([fn{ilp} '_feedforward']) = tempdata.(fn{ilp});
            data.([fn{ilp} '_feedback'])    = data.(fn{ilp});
        end
        data = rmfield(data,fn);
        return;
end
%-- channel selection
channel = ft_channelselection(channel, data.label);
chanreject = 1:nchn;
chanreject(match_str(data.label, channel)) = [];
chslidx = zeros(nchn,nchn);
chslidx(chanreject,:) = 1;   chslidx(:,chanreject) = 1;
chslidx = find(chslidx);                             % not tilde for nan
%-- combine
pickidx  =  unique([chslidx; pickidx]);

%-- pickup
strforpick = '(';
for idm =1:ndimdat;
    if idm~=1,  strforpick = [strforpick, ', '];    end
    if idm == dim_chancmb
        strforpick = [strforpick, 'pickidx'];
    else
        strforpick = [strforpick, ':'];
    end
end
strforpick = [strforpick, ')'];
for ilp =1:length(fn)
   if islogical(data.(fn{ilp}))
     eval(sprintf('data.(fn{ilp})%s = false;',strforpick));
   else
     eval(sprintf('data.(fn{ilp})%s = nan;',strforpick));
   end
end

%-- output
switch cmbrepresentation
    case 'sparse'
      nanparam = 1:ndimdat;  nanparam(dim_chancmb) = [];
      isnanidx = true(nchncmb,1);
      for ilp =1:length(fn)
         tmpparam = reshape(permute(data.(fn{ilp}),[dim_chancmb nanparam]),nchncmb,[]);
         isnanidx = and(isnanidx, min(isnan(tmpparam),[],2));
      end
      %-- remove all-NaN channel
      strforpick = '(';
      for idm =1:ndimdat;
          if idm~=1,  strforpick = [strforpick, ', '];    end
          if idm == dim_chancmb
              strforpick = [strforpick, 'isnanidx'];
          else
              strforpick = [strforpick, ':'];
          end
      end
      strforpick = [strforpick, ')'];
      for ilp =1:length(fn)
          eval(sprintf('data.(fn{ilp})%s = [];',strforpick));
      end
      data.labelcmb(isnanidx,:) = [];
        
    case 'full'
      data     = rmfield(data,'labelcmb');
      dimtok   = [dimtok(1:(dim_chancmb-1)) {'chan' 'chan'} dimtok((dim_chancmb+1):end)];
      sizparam = [sizparam(1:(dim_chancmb-1)), nchn, nchn, sizparam((dim_chancmb+1):end)];
      for ilp =1:length(fn)
         data.(fn{ilp}) = reshape(data.(fn{ilp}),sizparam);
      end
end
data.dimord     = sprintf('%s_',dimtok{:});  data.dimord(end)=[];
data.validlabel = true(size(data.label));
data.validlabel(chanreject) = false;

%-- not keep channel
if strcmp(keepchannel,'no')
    switch cmbrepresentation
        case 'full'
            ndimdat     = length(dimtok);
            strforpick = '(';
            for idm =1:ndimdat;
                if idm~=1,  strforpick = [strforpick, ', '];    end
                if strcmp(dimtok{idm},'chan')
                    strforpick = [strforpick, 'data.validlabel'];
                else
                    strforpick = [strforpick, ':'];
                end
            end
            strforpick = [strforpick, ')'];
            for ilp =1:length(fn)
                eval(sprintf('data.(fn{ilp}) = data.(fn{ilp})%s;',strforpick));
            end
    end
    data.label(~data.validlabel)      = [];
    data.validlabel(~data.validlabel) = [];
end
    

    






