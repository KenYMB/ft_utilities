function [A,Ve] = bs_ARmodel_multi(dat,order,winlen,winshift)
% BS_ARMODEL_MULTI Moving window for multivariate model
%
% Usage:
%   [A,Ve] = bs_ARmodel_multi(dat,order,winlen,winshift);
% 
%   dat: input file in matlab format;
%   order: model order
%   winlen: window length
%   winshift: shift length of window (default=1)
% 
% Output:  
%   A is the name of AR coefficient file
%	Ve is the name of AR noise file
% 
% See also: bs_ARmodel_bi, mov_mul_model, mov_bi_model, moving_window.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 18-Sep-2007 15:43:37$
% SHIS UT-Houston, Houston, TX 77030, USA.

% Using: opssmov (modified by Yuasa @20170322)

% 20170323 Yuasa: modified for parallel computation & enable winshift
% 20170403 Yuasa: name change 'bs_mov_mul_model' -> 'bs_ARmodel_multi'

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
%-- processing
[pts, chan, trl] = size(dat);

%-- toi
startp  = 1;    endp  = pts;
pts     = endp - startp + 1;
dat=dat(startp:endp,:,:);

 % save channel chan -ascii;
 % save trail trl -ascii;
 % save points pts -ascii;
 % save window winlen -ascii;
 % save order order -ascii;
 % save winshift winshift -ascii;
writedat(['dataset_' suffix '.bin'],dat);

%-- reduce disk I/O
if ispc
    status = eval(['dos(''' sprintf('opssmov dataset_%s.bin  A_%s Ve_%s %d %d %d %d %d %d',suffix,suffix,suffix,chan,trl,pts,winlen,order,winshift) ''');']);
else
    status = eval(['unix(''' sprintf('./opssmov dataset_%s.bin  A_%s Ve_%s %d %d %d %d %d %d',suffix,suffix,suffix,chan,trl,pts,winlen,order,winshift) ''');']);
end%if

%-- if error
if status ~= 0 || ~exist(['A_' suffix],'file') || ~exist(['Ve_' suffix],'file')
    error('Cannot compute AR-model with ''opssmov''!');
end%if

A=load(['A_' suffix]);
Ve=load(['Ve_' suffix]);

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
    delete(['A_' suffix])
    delete(['Ve_' suffix])

    %-- restore directory
    cd(cdir);
    warning(warnState);
end%bs_finish

end%bs_mov_mul_model
% [EOF]