function [dat,beta,x] = ft_pow_polyremoval(dat, order, begsample, endsample, flag, exinterp)

% FT_POW_POLYREMOVAL removed an Nth order polynomal from the data
% for power data.
%
% Use as
%   dat = ft_pow_polyremoval(dat, order, begsample, endsample, flag, exinterp)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      the order of the polynomial
%   begsample  index of the begin sample for the estimate of the polynomial
%   endsample  index of the end sample for the estimate of the polynomial
%   flag       apply polyremoval after separate the dat in two parts
%   exinterp   intepolate extravalues before applying polyremoval
%
% If begsample and endsample are not specified, it will use the whole
% window to estimate the polynomial.
%
% For example
%   ft_preproc_polyremoval(dat, 'spline')
% removes the  linear trend.
%
% See also FT_PREPROC_BASELINECORRECT, FT_PREPROC_DETREND

% Copyright (C) 2008-2014, Robert Oostenveld
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
% $Id: ft_preproc_polyremoval.m$

% Using: findpeaks, interp1

% 20170420 Yuasa

% take the whole segment if begsample and endsample are not specified
if nargin<3 || isempty(begsample)
  begsample = 1;
end
if nargin<4 || isempty(endsample)
  endsample = size(dat,2);
end
if nargin<5 || isempty(flag)
  flag = 0;
end
if nargin<6 || isempty(exinterp)
  exinterp = 0;
end

% this does not work on integer data
typ = class(dat);
dat = cast(dat, 'double');

% get min-peaks
[nchan, nsamples] = size(dat);

if exinterp
    padsamples  = ceil(nsamples*0.05);      % 5%
    if begsample ~= 1, begsample = begsample + padsamples; end
    if endsample ~= nsamples,
        endsample = endsample + padsamples;
    else
        endsample = endsample + padsamples*2;
    end
else
    padsamples = 0;
end
calcsamples = nsamples+padsamples*2;
basedat = zeros(nchan, calcsamples);

for j = 1:nchan
     tempdat = dat(j,:);
    if exinterp && length(tempdat)>=4 && sum(~isnan(tempdat(:)))>=4
        tempdat = interp1(1:nsamples,tempdat,(1:(nsamples+padsamples*2))-padsamples,'linear','extrap');
    else
        tempdat = interp1(1:nsamples,tempdat,(1:(nsamples+padsamples*2))-padsamples,'linear','extrap');
    end
    
    [minpeak, peaktp] = findpeaks(tempdat.*-1);
    minpeak = minpeak * -1;
    if numel(peaktp)<4      % get minimum for few min-peak
        mintp   = find(tempdat == nanmin(tempdat));
        nexdat  = tempdat;  nexdat(mintp) = nan;
        nmintp  = find(nexdat(:) == nanmin(nexdat(:)));
        peaktp  = unique([peaktp, mintp, nmintp]);
        minpeak = tempdat(peaktp);
    end
    if ~isempty(peaktp)          % example: NaN matrix
        if min(diff(tempdat(1:peaktp(1)))<=0)         % increase for timepoint 1
            minpeak = [tempdat(1:peaktp(1)), minpeak(2:end)];
            peaktp  = [1:peaktp(1), peaktp(2:end)];
        elseif peaktp(1)~=1                         % set first value
            minpeak = [tempdat(1), minpeak(1:end)];
            peaktp  = [1, peaktp(1:end)];
        end
        if min(diff(tempdat(peaktp(end):end))>=0)     % increase for timepoint end
            minpeak = [minpeak(1:(end-1)), tempdat(peaktp(end):end)];
            peaktp  = [peaktp(1:(end-1)), peaktp(end):calcsamples];
        elseif peaktp(end)~=calcsamples                % set last value
            minpeak = [minpeak(1:end), tempdat(end)];
            peaktp  = [peaktp(1:end), calcsamples];
        end
        while numel(peaktp)<4 % for interp1
            nexdat(nmintp) = nan;
            nmintp  = find(nexdat(:) == nanmin(nexdat(:)));
            peaktp  = unique([peaktp, nmintp]);
            minpeak = tempdat(peaktp);
        end
        basedat(j,:) = interp1(peaktp, minpeak, 1:calcsamples, 'pchip');
    else
        basedat(j,:) = tempdat;
    end
end

% construct a "time" axis
basis    = (1:calcsamples)-1;

% create a set of basis functions that will be fitted to the data (orders x timepoints)
if ~flag
  x = zeros(order+1,calcsamples);
  for i = 0:order
    x(i+1,:) = basis.^(i);
  end
else        % fit using a part of data
  x = zeros(order*(calcsamples-2)+1,calcsamples);
  x(1,:) =  basis.^0;
  for i = 1:order
    for j=1:(calcsamples-2)
      x((i-1)*(calcsamples-2)+j+1,:) = [zeros(1,calcsamples-j-2) basis(1:(j+2)).^(i)];
    end
  end
end

% estimate the contribution of the basis functions (nchan x orders)
beta    = basedat(:,begsample:endsample)*x(:,begsample:endsample)'/(x(:,begsample:endsample)*x(:,begsample:endsample)');

% remove the estimated basis functions
dat = dat - beta*x(:,(1:nsamples)+padsamples);

% convert back to the original input type
dat = cast(dat, typ);
