function ft_mnefif_write_coord_trans(output, trans, unit, input)
%
% FT_MNEFIF_WRITE_COORD_TRANS(output, trans [,unit] [,input])
%
% Writes a coordinate transformation structure
% if 'unit' is indicated, 'trans' is converted to meter-base value
% (MNEtoolbox interpret transformation information as "m" data)
% if 'input' is not indicated, prepare 'mne-trans.fif' in the same
% directory
%
% trans             % The coordinate transfomation structure
% trans.from        % Source coordinate system
% trans.to          % Destination coordinate system
% trans.trans       % The forward transform [4x4]
% 
% output   = 'filename';    % output file name
%
% -- option --
% unit     = 'mm'/'cm'/'m'(default); % unit of transform
% input    = 'filename';    % loaded file name (output file is created based on this file)
% 
% -- example --
% transMEG2MRI          = [];
% transMEG2MRI.from     = 4;  %FIFFV_COORD_HEAD
% transMEG2MRI.to       = 5;  %FIFFV_COORD_MRI
% transMEG2MRI.trans    = transform_vox2spm/transform_vox2neuromag;
% ft_mnefif_write_coord_trans(transform_filename, transMEG2MRI, 'mm');

% using: built-in

% 20160929 Yuasa: crate based on fiff_write_coord_trans
% 20161212 Yuasa: minor fix
% 20170814 Yuasa: change name from mnefif to neurofif

narginchk(2,inf)

if nargin==2
    ft_neurofif_write_coord_trans(output, trans);
elseif nargin==3
    ft_neurofif_write_coord_trans(output, trans, unit);
else
    ft_neurofif_write_coord_trans(output, trans, unit, input);
end

%{
me='original:ft_mnefif_write_coord_trans';

narginchk(2,4)

FIFF_COORD_TRANS=222;
FIFFT_COORD_TRANS_STRUCT=35;
FIFFV_NEXT_SEQ=0;

FIFF_FILE_ID=100;
FIFF_DIR_POINTER=101;
FIFF_FREE_LIST=106;
FIFF_NOP=108;

%?typedef struct _fiffCoordTransRec {
%  fiff_int_t   from;                   /*!< Source coordinate system. */
%  fiff_int_t   to;                     /*!< Destination coordinate system. */
%  fiff_float_t rot[3][3];              /*!< The forward transform (rotation part) */
%  fiff_float_t move[3];                /*!< The forward transform (translation part) */
%  fiff_float_t invrot[3][3];           /*!< The inverse transform (rotation part) */
%  fiff_float_t invmove[3];             /*!< The inverse transform (translation part) */
%} *fiffCoordTrans, fiffCoordTransRec;  /*!< Coordinate transformation descriptor */

if nargin < 3 || isempty(unit)
    unit = 'm';
end
if nargin < 4 || isempty(input)
    input = [fileparts(mfilename('fullpath')) '/mne-trans.fif'];
end
if length(output)<4 || ~(strcmp(output(end-3:end),'.fif') || strcmp(output(end-4:end),'.fiff'))
    output = [output '.fif'];
end
assert(logical(exist(input,'file')),me,'There is no template file');
if logical(exist(output,'file'))
    warning('The output file already exist! Now copying...');
    copyfile(output, [output '.bk']);
end

fid = fopen(input,'r');
sampledata = fread(fid,'int32','b');

%-- convert unit
if strcmp(unit,'mm')
    convuni = 10^-3;
elseif strcmp(unit,'cm')
    convuni = 10^-2;
else
    convuni = 10^0;
end

%-- get all tag.kind
tagidx = 1;  tagkind = [];
while tagidx < length(sampledata)
    tagkind = [tagkind; [sampledata(tagidx) tagidx]];
    if sampledata(tagidx+3)>0
        tagidx = sampledata(tagidx+3)/4 + 1;
    else
        tagidx = (tagidx+4) + sampledata(tagidx+2)/4;
    end
end

fclose(fid);
assert(tagkind(1,1)==FIFF_FILE_ID && tagkind(2,1)==FIFF_DIR_POINTER,me,'Input file is wrong');

startidx = tagkind(tagkind(:,1)==FIFF_COORD_TRANS,2);
endidx   = tagkind(tagkind(:,1)==FIFF_NOP,2);

%-- write pre-trans information
fid = fopen(output,'w+');

for idx=1:startidx-1
    count = fwrite(fid,int32(sampledata(idx)),'int32','b');
    if count ~= 1
        fclose(fid);
        error(me,'write failed');
    end
end    

%-- fiff write
datasize=4*2*12 + 4*2;
count = fwrite(fid,int32(FIFF_COORD_TRANS),'int32','b');
if count ~= 1
    fclose(fid);
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_COORD_TRANS_STRUCT),'int32','b');
if count ~= 1
    fclose(fid);
    error(me,'write failed');
end
count = fwrite(fid,int32(datasize),'int32','b');
if count ~= 1
    fclose(fid);
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFV_NEXT_SEQ),'int32','b');
if count ~= 1
    fclose(fid);
    error(me,'write failed');
end
%
%   Start writing fiffCoordTransRec
%
count = fwrite(fid,int32(trans.from),'int32','b');
if count ~= 1
    fclose(fid);
    error(me,'write failed');
end
count = fwrite(fid,int32(trans.to),'int32','b');
if count ~= 1
    fclose(fid);
    error(me,'write failed');
end
%
%   The transform...
%
rot=trans.trans(1:3,1:3)';
move=trans.trans(1:3,4)' .* convuni;
count = fwrite(fid,single(rot),'single','b');
if count ~= 9
    fclose(fid);
    error(me,'write failed');
end
count = fwrite(fid,single(move),'single','b');
if count ~= 3
    fclose(fid);
    error(me,'write failed');
end
%
%   ...and its inverse
%
trans_inv=inv(trans.trans);
rot=trans_inv(1:3,1:3)';
move=trans_inv(1:3,4)' .* convuni;
count = fwrite(fid,single(rot),'single','b');
if count ~= 9
    fclose(fid);
    error(me,'write failed');
end
count = fwrite(fid,single(move),'single','b');
if count ~= 3
    fclose(fid);
    error(me,'write failed');
end

%-- write pre-trans information
for idx=endidx:length(sampledata)
    count = fwrite(fid,int32(sampledata(idx)),'int32','b');
    if count ~= 1
        fclose(fid);
        error(me,'write failed');
    end
end

fclose(fid);
fprintf('Success output tansform-data\n');
%}
