function BS_ARMODEL_BI(dat,order,winlen,winshift,savedir)
% MOV_BI_MODEL Moving window for bivariate models
%
% Usage:
%   bs_ARmodel_bi(dat,order,winlen,winshift,savedir);
%
% Input(s):
%   dat      -   input file in matlab format;
%   order    -   model order
%   winlen   -   window length
%   winshift -   shift length of window (default=1)
%   savedir  -   directory name to save files (default='./Movingwindow_Coefficient'])
%
% Output(s):
%   in the bsmartroot/Movingwindow_Coefficient directory
% 
% See also: bs_ARmodel_multi, mov_bi_model, mov_mul_model, moving_window_pairwise.

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 09:51:48$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Using: opssmov (modified by Yuasa @20170322)

% 20170323 Yuasa: modified for parallel computation & enable winshift
% 20170403 Yuasa: name change 'bs_mov_bi_model' -> 'bs_ARmodel_bi'

narginchk(3,5);
%-- defaultvalue
if nargin<4 || isempty(winshift),   winshift = 1;    end
if nargin<5,                        savedir = '';    end

if isempty(savedir),    savedir='./Movingwindow_Coefficient';   end
[st, pathinfo] = fileattrib(savedir);
if st,  PathName = pathinfo.Name;
else
    mkdir(savedir);
    [st, pathinfo] = fileattrib(savedir);
    PathName = pathinfo.Name;
end

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
[pts, chanall, trl] = size(dat);
chan = 2;

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

[st, coefile] = fileattrib([PathName, filesep,'*']);
if st
  warning('delete all coefficients files in the directory');
  delete(coefile);
end

for i=1:(chanall-1)
    for j=(i+1):chanall
        dat1=dat(:,i,:);
        dat2=dat(:,j,:);
        dat3=cat(2,dat1,dat2);
        writedat(['dataset_' suffix '.bin'],dat3);
        %-- reduce disk I/O
        if ispc
            status = eval(['dos(''' sprintf('opssmov dataset_%s.bin  A_%s Ve_%s %d %d %d %d %d %d',suffix,suffix,suffix,chan,trl,pts,winlen,order,winshift) ''');']);
        else
            status = eval(['unix(''' sprintf('./opssmov dataset_%s.bin  A_%s Ve_%s %d %d %d %d %d %d',suffix,suffix,suffix,chan,trl,pts,winlen,order,winshift) ''');']);
        end%if
        
        %-- if error
        if status ~= 0
            error('Cannot compute AR-model between %d & %d channels with ''opssmov''!',i,j);
        end%if
        
        ii=num2str(i);jj=num2str(j);
        FileName=['AR_C_' ii 'and' jj];
        fA=fullfile(PathName,FileName);
        movefile(['A_' suffix],fA);
        FileName2=['AR_N_' ii 'and' jj];
        fVe=fullfile(PathName,FileName2);
        movefile(['Ve_' suffix],fVe);

        delete(['dataset_' suffix '.bin'])
    end
end

fprintf(1,'All model coefficient save in %s directory\n',PathName);

bs_finish;  % end of this function

catch ME
%-- clean job zero
delete(['dataset_' suffix '.bin'])
delete(['A_' suffix])
delete(['Ve_' suffix])

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

    %-- restore directory
    cd(cdir);
    warning(warnState);
end%bs_finish

end%bs_mov_bi_model
% [EOF]
