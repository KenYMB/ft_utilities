function AIC = bs_aic_test(dat,maxorder,winlen,winshift)
% AIC_TEST Compute the AIC
% 
% Usage:
%   [AIC] = bs_aic_test(data,maxorder,winlen,winshift)
% 
%   dat: data set in Matlab format
%   maxorder: the maximum model order you want to try.  (default=10)
%   winlen: window length 
%   winshift: shift length of window (default=1)
% 
%   AIC: each row is the AIC coefficient for each window.
% 
% See also: aic_test, aic.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 18-Sep-2007 15:41:17$
% SHIS UT-Houston, Houston, TX 77030, USA.

% Using: opssaic (modified by Yuasa @20170322)

% 20170322 Yuasa: modified for parallel computation & enable winshift

narginchk(2,4);
%-- defaultvalue
if nargin<3 || isempty(maxorder),   maxorder = 10;    end
if nargin<4 || isempty(winshift),   winshift = 1;    end

%-- parse directory
cdir = pwd;                 % find current directory
 % p = mfilename('fullpath');
p = which('bsmart');
fdir = fileparts(p);        % find function directory
cd(fdir);                   % change current dir to function dir
suffix  = sprintf('%d',fix(10*rand(1,6)));

warnState = warning('off','MATLAB:DELETE:FileNotFound');

try
%-- save opssaic input arguments
[pts, chan, trl] = size(dat);
 % save channel chan -ascii;
 % save trail trl -ascii;
 % save points pts -ascii;
 % save window winlen -ascii;
 % save order maxorder -ascii;
 % save winshift winshift -ascii;
writedat(['dataset_' suffix '.bin'],dat);

%-- reduce disk I/O
if ispc
    status = eval(['dos(''' sprintf('opssaic dataset_%s.bin  A_%s Ve_%s AIC_now_%s %d %d %d %d %d %d',suffix,suffix,suffix,suffix,chan,trl,pts,winlen,maxorder,winshift) ''')']);
else
    status = eval(['unix(''' sprintf('./opssaic dataset_%s.bin  A_%s Ve_%s AIC_now_%s %d %d %d %d %d %d',suffix,suffix,suffix,suffix,chan,trl,pts,winlen,maxorder,winshift) ''')']);
end%if

%-- if error
if status ~= 0 || ~exist(['AIC_now_' suffix],'file')
    error('Cannot compute AIC with ''opssaic''!');
end%if
%-- otherwise
AIC = load(['AIC_now_' suffix]);

bs_finish;  % end of this function

catch ME
bs_finish;  % end of this function
rethrow(ME);
end

function bs_finish
    %-- clean job one
     % delete channel
     % delete trail
     % delete points
     % delete window
     % delete order
     % delete winshift
    delete(['dataset_' suffix '.bin'])
    %-- clean job two
    delete(['AIC_now_' suffix])
    delete(['A_' suffix])
    delete(['Ve_' suffix])

    %-- restore directory
    cd(cdir);
    warning(warnState);
end%bs_finish

end%bs_aic_test
% [EOF]