function LE = bs_stability_test(arcoeff,arnoise,Tsamp,iswaitbar)
% Compute Lyapunov exponent
% 
% Usage: 
%   LE = bs_stability_test(arcoeff,arnoise,Tsamp);
%   LE = bs_stability_test(arcoeff,arnoise,Tsamp,'nowaitbar');
% 
%   Tsamp :   Number of samples to generate;
% 
% Note: 
%   In fact we do not use ARNOISE in Lyapunov calculation, 
%   but need it in MAR_make routine
%
% See also: lyap_batch, lyap, Lyapunov.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 10:39:05$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Hualou Liang, 02/09/99, FAU
% 

% Using: lyap

% 20170324 Yuasa: modified "lyap_batch" for 'winshift' using data <- no need to change?
% 20170405 Yuasa: add option to suppress waitbar

if nargin>=4 && ischar(iswaitbar) && strcmp(iswaitbar,'nowaitbar')
    iswaitbar = false;
else
    iswaitbar = true;
end

LE = [];
len = size(arnoise,1);
if iswaitbar
  hw = waitbar(0,'Please wait...','Name','Stability test','CreateCancelBtn','closereq','Visible','off');
  set(hw.Children(1),'String','OK','Visible','off');       set(hw,'CloseRequestFcn','','Visible','on');
else
  fprintf(1,'Stability Testing...    \n');
end
try
  for t = 1:len,
      if iswaitbar,     waitbar((t-1)/len);
      else              fprintf('\b\b\b\b\b% 3d%%\n',round((t-1)/len*100));
      end
      mar = MAR_make(arcoeff(t,:), arnoise(t,:));
      le = lyap(mar,Tsamp,1000);          % LE = lyap(mar, T, discardnum);
      LE = cat(2,LE,le);
  end;
  LE = LE';
  %-- finish waitbar
  if iswaitbar
    waitbar(1,hw,'Finish');    set(hw.Children(1),'Visible','on');
    set(hw,'CloseRequestFcn','closereq','HandleVisibility','on');
  end
  
catch ME
  %-- finish waitbar
  if iswaitbar
    waitbar(t/len,hw,'Error');    set(hw.Children(1),'Visible','on');
    set(hw,'CloseRequestFcn','closereq','HandleVisibility','on');
  else
    fprintf('\b\b\b\b\b% 3d%%\n',round(t/len*100));
  end
  
  rethrow(ME);
end
end%bs_stability_test

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
% [EOF]
