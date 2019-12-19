function [resid] = bs_whiteness_test(dat,order,winlen,winshift)
% WHITENESS_TEST Whiteness test
% 
% Usage:
%   [resid] = bs_whiteness_test(data,order,winlen,winshift)
% 
%   dat: data set in Matlab format
%   order: model order
%   winlen: window length
%   winshift: shift length of window (default=1)
%
%   resid: residuals probabilities
% 
% See also: whiteness_test, whiteness.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 22:43:55$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Using: opsswhite (modified by Yuasa @20170322)

% 20170322 Yuasa: modified for parallel computation & enable winshift

narginchk(3,4);
%-- defaultvalue
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
[pts, chan, trl] = size(dat);
 % save channel chan -ascii;
 % save trail trl -ascii;
 % save points pts -ascii;
 % save window winlen -ascii;
 % save order order -ascii;
 % save winshift winshift -ascii;
writedat(['dataset_' suffix '.bin'],dat);

%-- reduce disk I/O
if ispc
    status = eval(['dos(''' sprintf('opsswhite dataset_%s.bin  A_%s Ve_%s %d %d %d %d %d %d resid_%s.out',suffix,suffix,suffix,chan,trl,pts,winlen,order,winshift,suffix) ''')']);
else
    status = eval(['unix(''' sprintf('./opsswhite dataset_%s.bin  A_%s Ve_%s %d %d %d %d %d %d resid_%s.out',suffix,suffix,suffix,chan,trl,pts,winlen,order,winshift,suffix) ''')']);
end%if

%-- if error
if status ~= 0 || ~exist(['resid_' suffix '.out'],'file')
    error('Cannot test whiteness with ''opsswhite''!');
end%if
%-- otherwise
resid = load(['resid_' suffix '.out']);

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
    delete(['resid_' suffix '.out'])
    delete(['A_' suffix])
    delete(['Ve_' suffix])

    %-- restore directory
    cd(cdir);
    warning(warnState);
end%bs_finish

end%bs_whiteness_test
% [EOF]