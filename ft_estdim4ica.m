function [Dim] = ft_estdim4ica(cfg, data)

% FT_ESTDIM4ICA estimates ideal demenstions of pca for ica - especially for 'runica'.
% The estimation method is adopted in 'fastica'.
% 
% The rank is determined from the eigenvalues - and not directly by
% using the function rank - because function rank uses svd, which
% in some cases gives a higher dimensionality than what can be used
% with eig later on (eig then gives negative eigenvalues).
%
% Use as
%   [dim] = ft_estdim4ica(cfg, data)
%
% where the data comes from FT_PREPROCESSING and the configuration
% structure can contain
%   cfg.channel      = cell-array with channel selection (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.demean       = 'no' or 'yes', whether to demean the input data (default = 'yes')
%   cfg.doscale      = 'no' or 'yes', whether to scale the input data to approximately unity (default = 'yes')
% 
% Following is a sample usage of this function:
%     cfg                 = [];
%     cfg.updatesens      = 'yes';
%     cfg.demean          = 'yes';
%     cfg.doscale         = 'yes';
%     cfg.method          = 'runica';
%     cfg.runica.pca      = ft_estdim4ica(cfg,data);
%     ica      = ft_componentanalysis(cfg,data);
% 
% See also FT_COMPONENTANALYSIS

% using: fieldtrip

% 20160805 Yuasa: customize for ica: add components button


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on ft_componentanalysis & fastica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-- do the general setup of the function
ft_defaults;

%-- check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

%-- check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'blc', 'demean'});

%-- set the defaults
cfg.demean          = ft_getopt(cfg, 'demean',       'yes');
cfg.trials          = ft_getopt(cfg, 'trials',       'all', 1);
cfg.channel         = ft_getopt(cfg, 'channel',      'all');
cfg.doscale         = ft_getopt(cfg, 'doscale',      'yes');

%-- apply cfg
tmpcfg = [];
tmpcfg.trials = cfg.trials;
tmpcfg.channel = cfg.channel;
data = ft_selectdata(tmpcfg, data);

%-- get parameters
Ntrials  = length(data.trial);
Nchans   = length(data.label);
if Nchans==0
    error('no channels were selected');
end

Nsamples = zeros(1,Ntrials);
for trial=1:Ntrials
    Nsamples(trial) = size(data.trial{trial},2);
end

%-- apply demean
if strcmp(cfg.demean, 'yes')
    fprintf('baseline correcting data \n');
    for trial=1:Ntrials
        data.trial{trial} = ft_preproc_baselinecorrect(data.trial{trial});
    end
end

%-- apply doscale
if strcmp(cfg.doscale, 'yes')
    scale = norm((data.trial{1}*data.trial{1}')./size(data.trial{1},2));
    scale = sqrt(scale);
    if scale ~= 0
        fprintf('scaling data with 1 over %f\n', scale);
        for trial=1:Ntrials
            data.trial{trial} = data.trial{trial} ./ scale;
        end
    else
        fprintf('no scaling applied, since factor is 0\n');
    end
else
    fprintf('no scaling applied to the data\n');
end

%-- concatenating data
fprintf('concatenating data');
mixedsig = zeros(Nchans, sum(Nsamples));
for trial=1:Ntrials
    fprintf('.');
    begsample = sum(Nsamples(1:(trial-1))) + 1;
    endsample = sum(Nsamples(1:trial));
    mixedsig(:,begsample:endsample) = data.trial{trial};
end
fprintf('\n');
fprintf('concatenated data matrix size %dx%d\n', size(mixedsig,1), size(mixedsig,2));

%-- remean
mixedsig = mixedsig - mean(mixedsig,2) * ones (1,size(mixedsig, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on pcamat in fastica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-- default parameters
b_verbose     = true;
oldDimension  = size (mixedsig, 1);
firstEig      = 1;
lastEig       = oldDimension;
rankTolerance = 1e-7;

%-- Calculate the eigenvalues and eigenvectors of covariance matrix.
covarianceMatrix = cov(mixedsig', 1);
[E, D] = eig (covarianceMatrix);
eigenvalues = sort(diag(D),'descend');

%-- The rank is determined from the eigenvalues (to avoid eig giving negative eigenvalues)
maxLastEig = sum (diag (D) > rankTolerance);
assert(maxLastEig>0, ['Eigenvalues of the covariance matrix are' ...
                      ' all smaller than tolerance [ %g ].\n' ...
                      'Please make sure that your data matrix contains' ...
                      ' nonzero values.\nIf the values are very small,' ...
                      ' try rescaling the data matrix.\n'], rankTolerance);
if lastEig > maxLastEig,    lastEig = maxLastEig;   end

%-- Drop the smaller eigenvalues
if lastEig < oldDimension
    lowerLimitValue = (eigenvalues(lastEig) + eigenvalues(lastEig + 1)) / 2;
else
    lowerLimitValue = eigenvalues(oldDimension) - 1;
end
lowerColumns = diag(D) > lowerLimitValue;

%-- Drop the larger eigenvalues
if firstEig > 1
    higherLimitValue = (eigenvalues(firstEig - 1) + eigenvalues(firstEig)) / 2;
else
    higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(D) < higherLimitValue;

%-- Combine the results from above
selectedColumns = lowerColumns & higherColumns;

Dim = sum(selectedColumns);

%-- print some info for the user
if b_verbose
  fprintf ('Selected [ %d ] dimensions.\n', sum (selectedColumns));
end
if sum (selectedColumns) ~= (lastEig - firstEig + 1),
  error ('Selected a wrong number of dimensions.');
end

if b_verbose
  fprintf ('Smallest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(lastEig));
  fprintf ('Largest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(firstEig));
  fprintf ('Sum of removed eigenvalues [ %g ]\n', sum(diag(D) .* ...
    (~selectedColumns)));
end
