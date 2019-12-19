function [s, cfg] = ft_statfun_bootstrap(cfg, dat, design)

% FT_STATFUN_BOOTSTRAP calculates the bootstrap test on the biological data in dat
% (the one-sample variable).
% Don't use with "cfg.correctm = 'cluster'"
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_bootstrap'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_depsamplesT(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (ivar) and the unit-of-observation (uvar) 
%          factor,  Nfac x Nreplications
%
% Configuration options
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification
%   cfg.ivar  = row number of the design that contains the labels of the conditions that must be 
%               compared (default=1). The labels are the numbers 1.

% Copyright (C) 2006, Eric Maris
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

% 20170428 Yuasa: only "one-sample bootstrap" is implemented

% set defaults
if ~isfield(cfg, 'computestat'),    cfg.computestat    = 'yes'; end 
if ~isfield(cfg, 'computecritval'), cfg.computecritval = 'no';  end
if ~isfield(cfg, 'computeprob'),    cfg.computeprob    = 'no';  end
if ~isfield(cfg, 'alpha'),          cfg.alpha          = 0.05;  end
if ~isfield(cfg, 'tail'),           cfg.tail           = 1;     end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') || strcmp(cfg.computecritval,'yes')
    error('This function can only output statistical value.');
end;

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
nreplc1 = sum(~isnan(dat(:,sel1)), 2);
nreplc2 = sum(~isnan(dat(:,sel2)), 2);
nrepl   = nreplc1 + nreplc2;
if any(nrepl<size(design,2)),
  ft_warning('Not all replications are used for the computation of the statistic.');
end;
%if nrepl<3
%    error('The data must contain at least three trials/subjects.');
%end;
[nsmpls, nunits] = size(dat);

if strcmp(cfg.computestat, 'yes')
  % compute the statistic
  if nargout ~= 1   % skip bootstrap for first step
      bootdat = dat;
  else
      bootidx = ceil(nsmpls*rand(nsmpls,nunits));
      bootidx = bootidx + repmat((0:(nunits-1))*nsmpls,nsmpls,1);
      bootdat = dat(bootidx);
  end
      
  if ~isempty(sel1)
    avg1 = nanmean(bootdat(:,sel1), 2);
  else
    avg1 = zeros(nsmpls,1);
  end
  if ~isempty(sel2)
    avg2 = nanmean(bootdat(:,sel2), 2);
  else
    avg2 = zeros(nsmpls,1);
  end
  
  s.stat = avg1 - avg2;
end

if strcmp(cfg.computecritval,'yes')
  % also compute the critical values
end

if strcmp(cfg.computeprob,'yes')
  % also compute the p-values
end

