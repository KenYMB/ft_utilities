function ratio = bs_consistency_test(dat,arcoeff,arnoise,winlen,winshift,fsample)
% CONSISTENCYTEST Consistency test
% 
% Sytax:
%   ratio = bs_consistency_test(data,A,Ve,winlen,winshift)
% 
%   data: data set in Matlab format
%   winlen: window length ([time point]/[s])
%   winshift: shift length of window ([time point]/[s])
%             if empty, estimated by data size, window length, and AR model size
%   fsample: sampling frequency
%            this value is used to decide lag of cross-correlation in "bs_constchk"
%            if this value is input, 'winlen', 'winshift' are interpreted
%            as second-order values.
%   -- following parameters must be calculated by above values --
%   A: AR coefficient (movedwindow x ((order+1)*chans^2))
%   Ve: AR noise      (movedwindow x (chans^2))
% 
% Ouput(s): 
%   ratio: consistency ratio in Matlab format
% 
% See also: consistencytest, movingwin, constchk, consistency_test.

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 15-Sep-2007 09:04:13$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Using: gtm_dist

% 20170323 Yuasa: modified "consistencytest" for 'winshift' using data
% 20170324 Yuasa: modify 'nlags'
% 20170405 Yuasa: bug fix: MAR_gen -> "./lag"

narginchk(4,6);

%-- processing
[pts, chan, trl] = size(dat);

if nargin == 6
    winlen   = round(winlen .* fsample);
    winshift = round(winshift .* fsample);
    
    nlags    = round(fsample ./ 40);   % 5 for Fs=200Hz
else
    nlags    = [];
end

%-- data check
if chan ~= sqrt(size(arnoise,2))
    error('indicated data and AR models are inconsistent');
end
if nargin<5 || isempty(winshift),
    winshift  = fix((pts - winlen)/size(arcoeff,1));
end

[R1,R2,diffR,ratio] = bs_movingwin(dat,arcoeff,arnoise,chan,trl,pts,winlen,winshift,nlags);

ratio = 100*(1-ratio);
ratio(ratio<0) = 0; % neg val correction (if coeff for AR-sim data >> input data)

end%bs_consistency_test

function [R1,R2,diffR,ratio]=bs_movingwin(gelb,arcoeff,arnoise,chan,trl,pts,winlen,winshift,nlags)
% Input:
%   gelb: T x chans x trials; T - length of the data
%   arcoeff: movedwindow x ((order+1)*chans^2)
%   arnoise: movedwindow x (chans^2)
%   nlags: calculate correlation between -nlags~nlags (default = 5)
% Output:
%   R1: radius of data1
%   R2: radius of data2
%   diffR: difference between R1 and R2
%   ratio: diffR/R1

if ~exist('nlags','var') || isempty(nlags),  nlags = 5;  end
nlags = max(1,abs(nlags));

shift=zeros(trl,1);
shiftmin=min(shift);  shfmax=max(shift);

maxnwinshift = size(arcoeff,1);

R1=[];
R2=[];
diffR=[];
ratio=[];
idx = 1;
for rec=0:winshift:(pts-winlen-shfmax)
  if idx <= maxnwinshift
%     dat1 = zeros(winlen, chan, trl);  % monkey data
%     for i=1:trl 
% 	    dat1(:,:,i) = gelb(rec+shift(i)+1:rec+shift(i)+winlen,:,i);
%     end
    dat1 = gelb(rec+shiftmin+1:rec+shiftmin+winlen,:,:);
    
    mar = MAR_make(arcoeff(idx,:), arnoise(idx,:));
    %-- generate simulation data
    dat2 = zeros(winlen, chan, trl);
    for itrl=1:trl,
        dat2(:,:, itrl) = MAR_gen(mar, winlen)';	% MAR_gen: chan x winlen
    end
    %-- check
    [R11, R21, diffR1, ratio1] = bs_constchk(dat1, dat2, nlags);
    
    R1=[R1;R11];
    R2=[R2;R21];
    diffR = [diffR;diffR1];
    ratio =[ratio;ratio1];
    
    idx = idx+1;
  else
    warnstat=warning('off','backtrace');
    warning('indicated data and calculated numbers of AR models are inconsistent');
    warning(warnstat);
    return;
  end
end
end%bs_movingwin

function [mar] = MAR_make(coeff,noise)
%
% To make it more general for any channel numbers
%
% A batch version, to covert MAR output to MAR structure
%
% Usage:
%    [mar] = MAR_make(coeff,noise);
% Example: 
%    [mar] = MAR_make(coeff,noise);
%    [x] = ensembding(mar, T, chans, N)
% Input: 
%    coeff: row vector (1 x ((order+1)*chans^2))
%               [O0R1C1,O0R1C2,...,O0R1Cn,O0R2C1,...,O1R1C1,...,OmRnCn]
%    noise: row vector (1 x (chans^2))
% MAR is modeled as 
%   A0*X_t + A1*X_(t-1) + ... + A5*X_(t-5) = Et 
%
% Output:
%    mar: structure
%     mar.order               order of model
%     mar.lag(k).a            coeff matrix at lag k
%     mar.noise_cov           noise covariance matrix
% 

%-- generate structure MAR
mar = [];
mar.order = length(coeff)/length(noise) - 1;
chan = sqrt(length(noise));

for iO=1:mar.order,
  mar.lag(iO).a =  reshape(coeff(iO*chan*chan+1:(iO+1)*chan*chan),chan,chan)'; 
end

%-- Noise coefficient matrix
mar.noise_cov = reshape(noise,chan,chan)';

end% MAR_make

function [Y] = MAR_gen(mar, T, discardnum)
%
%   [Y] = mar_gen (mar, T)
%  Generate T samples of time series from MAR model
%   Y_t + A1*Y_(t-1) + ... =  C*W_t
% 
% Input: 
%   mar.lag(k).a  : is ar coefficient matrix at lag k
%   mar.noise_cov : estimated noise covariance
%   T             : number of samples to generate
%   discardnum    : discard the transient data length
%
% Output:
%      Y:  dim x T              
%

if nargin < 3,
  discardnum = 10^3;
end

lag = mar.order;  % Order of AR model
dim = size(mar.lag(1).a,1);

randvec = gsamp(T+discardnum,zeros(1,dim),mar.noise_cov)';  % dim x T+discardnum
%-- Generate first lag elements
Y = randvec;

%-- Generate rest of series
for i=lag+1:T+discardnum,
  for k=1:lag,
    Y(:,i) = Y(:,i) - mar.lag(k).a*Y(:,i-k);
  end
end

Y = Y(:,discardnum+1:discardnum+T);
end% MAR_gen

function x = gsamp(nsamp, mu, covar)
%GSAMP	Sample from a Gaussian distribution.
%
%	Description
%
%	X = GSAMP(NSAMP, MU, COVAR) generates a sample of size NSAMP from a
%	D-dimensional Gaussian distribution. The Gaussian density has mean
%	vector MU and covariance matrix COVAR, and the matrix X has NSAMP
%	rows in which each row represents a D-dimensional sample vector.
%
%
%	Copyright (c) Hualou Liang

d = size(covar, 1);
mu = reshape(mu, 1, d);   % Ensure that mu is a row vector
[evec, eval] = eig(covar);
coeffs = randn(nsamp, d)*sqrt(eval);

x = ones(nsamp, 1)*mu + coeffs*evec';
end%gsamp

function [R1, R2, diffR, ratio] = bs_constchk(dat1, dat2, nlags)
% For each data window, computing all chan-pairwise cross correlation
% function at different time lag, then check the ratio of radii
% Usage:
%   [R1, R2, diffR, ratio] = bs_constchk(data1, data2, nlags)
%   
% Input:
%   data1: sample data
%   data2: smimulation data based on MAR, same format as data1
%   nlags: calculate correlation between -nlags~nlags
% Output:
%   R1: radius of data1
%   R2: radius of data2
%   diffR: difference between R1 and R2
%   ratio: diffR/R1
%
if nargin<3,  nlags = [];  end

[val1] = corrlag(dat1,nlags);
[val2] = corrlag(dat2,nlags);
  

%-- estimate distance from each row to origin(0)
R1 = sqrt(sum(val1.^2));  % 1 x 1
R2 = sqrt(sum(val2.^2));  % 1 x 1

diffR = sqrt(gtm_dist(val1, val2));  % 1 x 1

ratio = diffR/R1;

end %bs_constchk

function [val] = corrlag(dat,nlags)
% calculate cross-correlation for various lags, and put in one row vector
% Input:
%  dat  : timepoints x channels x trials
%  nlags: calculate correlation between -nlags~nlags (default = 5)
% Output:
%  val: 1 x Len
%

if nargin==1 || isempty(nlags),  nlags = 5;  end
nlags = max(1,abs(nlags));

[pts, chan, trl] = size(dat);   % trl is size of pool

sel = nonzeros(tril(reshape(1:chan*chan, chan, chan),-1)); % select & remove autocorrelation 

  val = 0;
  for j=1:trl,
	val = val + xcorr(dat(:,:,j), nlags,'coeff');  % (2*nlags+1) x chan^2
  end  
  val = val/trl;
  val = val(:,sel);
  val = val(:)';     % runs (1) x p [No_of_variables(NoV)]  
end %corrlag
% [EOF]