function [res] = ft_mnefif_compute_raw_inverse(data,fname_inv,nave,lambda2,dSPM,sLORETA,labelfile)
%
% [res] = FT_MNEFIF_COMPUTE_RAW_INVERSE(data,fname_inv,nave,lambda2,dSPM,sLORETA,labelfile)
%
% An example on how to compute a L2-norm inverse solution
% Actual code using these principles might be different because 
% the inverse operator is often reused across data sets.
%
%
% data         - structure including "trial, time, label, fsample"
% data.trial   - Raw data {trials}[chs x times]
% data.time    - Time points {trials}[times]  (should be [s])
% data.label   - Sensor names {chs}
% data.fsample - Sampling frequency
% fname_inv    - Inverse operator file name
% nave         - Number of averages (scales the noise covariance)
%                If negative, the number of averages in the data will be used
% lambda2      - The regularization factor
%                lambda^2 = snr^-2, where snr is an estimate for amplitude SNR
%                default value of snr = 3
% dSPM         - do dSPM?
% sLORETA      - do sLORETA?
% labelfile    - The name of the label file (option)
%                If multiple labels are input as cellstr, logical OR is applied
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2008/05/26 10:49:26  msh
%   Update to incorporate the already weighted lead field basis
%
%   Revision 1.1  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%

% using: fieldtrip(, ft_hastoolbox, mne, match_str)

% 20161212 Yuasa: create based on mne_ex_compute_inverse.m
% 20161227 Yuasa: help update
% 20170208 Yuasa: avoid error if there is no available vertex
% 20170214 Yuasa: convert sparse matrix of res.sol to full matrix, 
%                 even if only a vertex is included
% 20170216 Yuasa: enable multiple labels
% 20170810 Yuasa: minor update

curpath = path;

me='original:ft_mnefif_compute_raw_inverse';

ft_defaults;
ft_hastoolbox('mne',1,0);

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

tic;
try
    
narginchk(4,7);

if nargin < 5
   dSPM = false;
end
if nargin < 6
   sLORETA = false;
end
%
%   Check sampling frequency
%
try
  if ~isfield(data,'fsample')
   if iscell(data.time)
     data.fsample = (length(data.time{1})-1) ./ (data.time{1}(end)-data.time{1}(1));
   else
     data.fsample = (length(data.time)-1) ./ (data.time(end)-data.time(1));
   end
  end
catch
    error(me,'Sampling frequency cannot be estimated');
end
%
%   Read the inverse operator
%
inv = mne_read_inverse_operator(fname_inv);
%
%   Set up the inverse according to the parameters
%
if nave < 0
    nave = 1;
end
inv = mne_prepare_inverse_operator(inv,nave,lambda2,dSPM,sLORETA);
%
%   Read label information
%
vertno = [];
for ihemi = 1:length(inv.src)
    vertno = [vertno reshape(inv.src(ihemi).vertno,1,[])-1];    % -1
end
if nargin >= 7
  vertices  = [];
  if ~iscell(labelfile),  labelfile = {labelfile};    end
  for ilabel = 1:length(labelfile)
    label     = mne_read_label_file(labelfile{ilabel});
    vertices  = [vertices reshape(intersect(vertno,label.vertices),1,[])];
    res.label{ilabel,1}     = label;
    res.label{ilabel}.path  = labelfile{ilabel};
  end
    selvert   = max(repmat(vertno',1,length(vertices)) == repmat(vertices,length(vertno),1),[],2);  % judge intersect
    vertno    = vertno(selvert);
else
    selvert   = true(inv.nsource, 1);
end
seltris   = reshape(repmat(selvert,1,3)',[],1);     % convert to triplets
%
%   Pick the correct channels from the data
%
selchan = match_str(data.label,inv.noise_cov.names);
assert(length(selchan)==length(inv.noise_cov.names),me,'Indicated inverse operator is not for this data');
for itrl=1:length(data.trial)
    data.trial{itrl} = data.trial{itrl}(selchan,:);
end
fprintf(1,'Picked %d channels from the data\n',length(selchan));

sols = cell(1,length(data.trial));
%
%   Check the available vertices
%
if isempty(selvert)
    warning(me, 'There is no available vertex.');
    for itrl=1:length(data.trial)
        sols{itrl} = [];
    end
else
    for itrl=1:length(data.trial)
        fprintf(1,'Computing inverse...');
        %
        %   Simple matrix multiplication followed by combination of the 
        %   three current components
        %
        %   This does all the data transformations to compute the weights for the
        %   eigenleads
        %   
        trans = diag(sparse(inv.reginv))*inv.eigen_fields.data*inv.whitener*inv.proj*double(data.trial{itrl});
        %
        %   Transformation into current distributions by weighting the eigenleads
        %   with the weights computed above
        %
        if inv.eigen_leads_weighted
           %
           %     R^0.5 has been already factored in
           %
           fprintf(1,'(eigenleads already weighted)...');
           sol   = inv.eigen_leads.data(seltris,:)*trans;
        else
           %
           %     R^0.5 has to factored in
           %
           fprintf(1,'(eigenleads need to be weighted)...');
           sol   = diag(sparse(sqrt(inv.source_cov.data(seltris,:))))*inv.eigen_leads.data(seltris,:)*trans;
        end

        if inv.source_ori == FIFF.FIFFV_MNE_FREE_ORI
            fprintf(1,'combining the current components...');
            sol1 = zeros(size(sol,1)/3,size(sol,2));
            for k = 1:size(sol,2)
                sol1(:,k) = sqrt(mne_combine_xyz(sol(:,k)));
            end
            sol = sol1;
        end
        if dSPM
            fprintf(1,'(dSPM)...');
            sol = inv.noisenorm(selvert,selvert)*sol;
        elseif sLORETA
            fprintf(1,'(sLORETA)...');
            sol = inv.noisenorm(selvert,selvert)*sol;
        end
        sols{itrl} = full(sol);
        fprintf(1,'trial%03d finish\n',itrl);
    end
end

res.inv    = inv;
res.sol    = sols;
res.vertno = vertno;
if isfield(data,'time')
 res.time  = data.time;
  if iscell(data.time)
   res.tmin  = data.time{1}(1);
  else
   res.tmin  = data.time(1);
  end
end
res.tstep = 1/data.fsample;
fprintf(1,'[done]\nTotal time is %.0f seconds\n',toc);

path(curpath);
return;

catch ME
    path(curpath);
    rethrow(ME);
end
end