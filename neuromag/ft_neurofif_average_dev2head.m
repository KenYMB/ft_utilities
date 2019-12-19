function ft_neurofif_average_dev2head(inputs, varargin)

% FT_NEUROFIF_AVERAGE_DEV2HEAD(inputs [,output] [,'rot', rottype])
% 
% output            % name of output file
%                   % default name is 'dev2head_ave.fiff'
% rot               % output type of rotation matrix
%                   % 'eye'(default): output identity matrix
%                   % 'average'     : output average across input data
% 
% ---- unusual options --------
% avgtype           % take average on ...
%                   % 'head2dev'(default): average head2dev matrix
%                   % 'dev2head'         : average dev2head matrix
% 
% 
% FT_NEUROFIF_AVERAGE_DEV2HEAD(inputs) just outputs Identity matrix with
% averaged head position.
% 

% using: fieldtrip(, ft_hastoolbox, mne), ft_neurofif_write_coord_trans,
%        cellfind, rmempty, rotmat3D, invrotmat3D

% 20170814 Yuasa
% 20180306 Yuasa: update the computation for rottype='eye'
%                (cuz, origin of translation matrix based on output origin)
% 20180330 Yuasa: change average method
%                (taking average of rotation matrix destroy orthogonality)


me='original:ft_neurofif_average_dev2head';

narginchk(1,inf)

ft_defaults;
ft_hastoolbox('mne',1,0);

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

%-- set 'rottype'
isrotopt = cellfind(varargin,'rot');
rottype = [];
if ~isempty(isrotopt) && isrotopt(1) < length(varargin)
    rottype = varargin{isrotopt(1)+1};
end
if strcmpi(rottype,'average') || strcmpi(rottype,'avg') || strcmpi(rottype,'ave')
    rottype = true;             % average
else
    rottype = false;            % eye
end

%-- set 'avgtype'
isavgopt = cellfind(varargin,'avgtype');
avgtype = [];
if ~isempty(isavgopt) && isavgopt(1) < length(varargin)
    avgtype = varargin{isavgopt(1)+1};
end
if strcmpi(avgtype,'dev2head')
    avgtype = false;            % dev2head
else
    avgtype = true;             % head2dev
end

%-- set 'output'
isoutput = 1:length(varargin);
rmidx    = [isrotopt; isrotopt+1];
rmidx(rmidx > length(varargin)) = [];
isoutput(rmidx) = [];
if isempty(isoutput)
    output   = 'dev2head_ave.fiff';
else
    output   = varargin{isoutput(1)};
    if iscell(output),   output = output{1};  end
end

%-- make load list
if ~iscell(inputs)
    inputs  = {inputs};
end
inputsInt = {};
for ilp=1:length(inputs)
  if ispc
    dirname   = fileparts(inputs{ilp});
    tmplist   = ls(inputs{ilp});
    tmplist   = cellstr([repmat([dirname filesep],size(tmplist,1),1) tmplist]);
  else
    tmplist   = reshape(strsplit(ls(inputs{ilp})),[],1);
  end
  inputsInt = [inputsInt; tmplist];
end
inputsInt = rmempty(inputsInt);
if isempty(inputsInt)
    error(message('MATLAB:open:fileNotFound', name));
end

%-- load input
fprintf(1,'Now loading...');
clear trans_dev2head
trans_dev2head(1:length(inputsInt)) = struct('from',[],'to',[],'trans',[]);
validlist = true(1,length(inputsInt));
for ilp=1:length(inputsInt)
    %-- load trans
    input = inputsInt{ilp};    
    trans_dev2head(ilp) = fiff_read_coord_trans(input);
    
    %-- check trans
    if trans_dev2head(ilp).from~=FIFF.FIFFV_COORD_DEVICE || trans_dev2head(ilp).to~=FIFF.FIFFV_COORD_HEAD
        fprintf(2,'transform matrix in ''%s'' is invalid.\n',input);
        trans_dev2head(ilp) = [];
        validlist(ilp)      = false;
    end
end
%-- apply inverse & compute rotation property & get the position of head origin
valididx  = find(validlist);
nvalid    = length(valididx);
rotaxis  = zeros(nvalid,3);
rottheta = zeros(nvalid,1);
headpos  = zeros(nvalid,3);
for ilp=1:nvalid
    if avgtype,    trans4avg = inv(trans_dev2head(valididx(ilp)).trans);    % 'head2dev'
    else           trans4avg = trans_dev2head(valididx(ilp)).trans;         % 'dev2head'
    end
    [rotaxis(ilp,:), rottheta(ilp)] = invrotmat3D(trans4avg(1:3,1:3));
    if  rottheta(ilp)<0
        rotaxis(ilp,:) = -rotaxis(ilp,:);
        rottheta(ilp)  = -rottheta(ilp);
    end
    headpos(ilp,:) = trans4avg(1:3,4);
end

%-- main rotation axis
repaxis = sum(rotaxis.*repmat(abs(rottheta),1,3),1)./sum(abs(rottheta));
for ilp=1:nvalid
    if rotaxis(ilp,:)*repaxis' < 0
        rotaxis(ilp,:) = -rotaxis(ilp,:);
        rottheta(ilp)  = -rottheta(ilp);
    end
end

%-- average trans
avgtransmat = zeros(4);     avgtransmat(4,4) = 1;
avgtransmat(1:3,4) = mean(headpos,1);
avgtransmat(1:3,1:3) = rotmat3D(sum(rotaxis.*repmat(abs(rottheta),1,3),1), mean(rottheta));
if ~avgtype,    avgtransmat     = inv(avgtransmat);         end     % false = 'dev2head', need inverse
if ~rottype,    avgtransmat(1:3,1:3) = eye(3);              end     % false = 'eye'

%-- make output data
avg_trans = struct('from',FIFF.FIFFV_COORD_DEVICE,'to',FIFF.FIFFV_COORD_HEAD,'trans',zeros(4));
avg_trans.trans = inv(avgtransmat);

%-- write out
fprintf(1,'writng...\n');
ft_neurofif_write_coord_trans(output, avg_trans, 'm',inputsInt{find(validlist,1)});

function [trans_head2mri] = fiff_read_coord_trans(transfile)

global FIFF;

[fid,~,dir] = fiff_open(transfile);
a = find([dir.kind] == FIFF.FIFF_COORD_TRANS);
pos = dir(a(1)).pos;
tag = fiff_read_tag(fid,pos);
trans_head2mri = tag.data;
fclose(fid);                    % original function forget 'fclose'

