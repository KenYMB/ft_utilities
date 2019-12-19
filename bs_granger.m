function [Fxy,Fyx] = bs_granger(dat,order,winlen,winshift,fsample,freq,iswaitbar)
% MOV_BI_GA Compute the granger causality from the moving window Bivariate models
% 
% Usage:
%   [Fxy,Fyx] = bs_granger(data,order,winlen,winshift,fsample,freq) 
%   [Fxy,Fyx] = bs_granger(data,order,winlen,winshift,fsample,freq,'nowaitbar') 
% 
% Input(s):
%   data     - data set in Matlab format
%   order    - model order
%   winlen   - window length [s]
%   winshift - shift length of window [s]
%   fsample  - Sampling frequency
%   freq     - a vector of frequencies of interest, usually   freq = 0:fs/2;
% 
% Output(s):
%   Fxy     - the causality measure from x to y
%   Fyx     - the causality measure from y to x
%             The order of Fx2y/Fy2x is 1 to 2:N; 2 to 3:N; ...; N-1 to N,
%              where N is the number of channels.
%             That is,
%              1st row=1&2; 2nd=1&3; ...; (N-1)th=1&N; ...; (N(N-1)/2)th=(N-1)&N.
% 
% See also: mov_bi_ga, one_bi_ga, pwcausal, pairwise_granger_causality.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 15:53:00$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Using: pwcausal (,armorf, spectrum)

% 20170324 Yuasa: modified "mov_bi_ga" to enable winshift
% 20170328 Yuasa: fix output index
% 20170405 Yuasa: add option to suppress waitbar

narginchk(4,7);
if nargin>=7 && ischar(iswaitbar) && strcmp(iswaitbar,'nowaitbar')
    iswaitbar = false;
else
    iswaitbar = true;
end

%-- processing
[pts, chan, trl] = size(dat);

%-- convert to time point
winlen    = round(winlen .* fsample);
winshift  = round(winshift .* fsample);

%-- toi
startp  = 1;    endp  = pts;
pts     = endp - startp + 1;

nwinshift   = ceil((pts-winlen+1)/winshift);
nchancmb    = (chan-1)*chan/2;
total       = nwinshift * nchancmb;
Fxy     = zeros(nchancmb,length(freq),nwinshift);
Fyx     = zeros(nchancmb,length(freq),nwinshift);
count   = 0;
if iswaitbar
  hw = waitbar(0,'Please wait...','Name','Granger Causality Calculation','CreateCancelBtn','closereq','Visible','off');
  set(hw.Children(1),'String','OK','Visible','off');       set(hw,'CloseRequestFcn','','Visible','on');
else
  fprintf(1,'Calculating Granger Causality...    \n');
end
try
  for iwin = 1:winshift:pts-winlen+1
      kch = 1;      % output index
      if iswaitbar,     waitbar(count/total);
      else              fprintf('\b\b\b\b\b% 3d%%\n',round(count/total*100));
      end
      for ich = 1:(chan-1)
          for jch = (ich+1):chan
              %-- reshape & pickup data [chans_of_interest x (data_points_of_interest x trials)]
              tmpdat = reshape(permute(dat(iwin+startp-1:iwin+startp+winlen-2,[ich,jch],:),[2 1 3]),2,[]);
              %-- main
              [pp,cohe,Fxy(kch,:,ceil(iwin/winshift)),Fyx(kch,:,ceil(iwin/winshift))]...
                  = pwcausal(tmpdat,trl,winlen,order,fsample,freq);
              
              count = count + 1;
              kch   = kch + 1;
          end
      end
  end
  
  %-- finish waitbar
  if iswaitbar
    waitbar(1,hw,'Finish');    set(hw.Children(1),'Visible','on');
    set(hw,'CloseRequestFcn','closereq','HandleVisibility','on');
  else
    fprintf('\b\b\b\b\b% 3d%%\n',round(count/total*100));
  end
  
catch ME
  %-- finish waitbar
  if iswaitbar
    waitbar(count/total,hw,'Error');    set(hw.Children(1),'Visible','on');
    set(hw,'CloseRequestFcn','closereq','HandleVisibility','on');
  end
  
  rethrow(ME);
end

end%bs_granger
% [EOF]