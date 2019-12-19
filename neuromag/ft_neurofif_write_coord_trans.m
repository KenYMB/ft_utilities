function ft_neurofif_write_coord_trans(output, trans, unit, input)
%
% FT_NEUROFIF_WRITE_COORD_TRANS(output, trans [,unit] [,input])
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
% trans.trans       % The forward transform in fieldtrip manar [4x4]
% 
% output   = 'filename';    % output file name
%
% -- option --
% unit     = 'mm'/'cm'/'m';     % unit of transform matrix (to convert into 'm')
%                               % (-trans.fif is interpreted as 'm' unit)
% -- unusual option --
% input    = 'filename';        % loaded file name (fileID is copied from this file)
% 
% -- example --
% transMEG2MRI          = [];
% transMEG2MRI.from     = 4;  %FIFFV_COORD_HEAD
% transMEG2MRI.to       = 5;  %FIFFV_COORD_MRI
% transMEG2MRI.trans    = transform_vox2spm/transform_vox2neuromag;
% ft_neurofif_write_coord_trans(transform_filename, transMEG2MRI, 'mm');

% using: fieldtrip(, ft_hastoolbox, mne, ft_estimate_units)

% 20160929 Yuasa: crate based on fiff_write_coord_trans
% 20161212 Yuasa: minor fix
% 20170810 Yuasa: update not to require input file
% 20170814 Yuasa: change name from mnefiff to neurofif

me='original:ft_neurofif_write_coord_trans';

narginchk(2,4)

ft_defaults;
ft_hastoolbox('mne',1,0);

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

%?typedef struct _fiffCoordTransRec {
%  fiff_int_t   from;                   /*!< Source coordinate system. */
%  fiff_int_t   to;                     /*!< Destination coordinate system. */
%  fiff_float_t rot[3][3];              /*!< The forward transform (rotation part) */
%  fiff_float_t move[3];                /*!< The forward transform (translation part) */
%  fiff_float_t invrot[3][3];           /*!< The inverse transform (rotation part) */
%  fiff_float_t invmove[3];             /*!< The inverse transform (translation part) */
%} *fiffCoordTrans, fiffCoordTransRec;  /*!< Coordinate transformation descriptor */

if nargin < 3 || isempty(unit)
    unit = ft_estimate_units(norm(trans.trans(1:3,4)));     % estimate from move vector
end
if length(output)<5 || ~(strcmp(output(end-3:end),'.fif') || strcmp(output(end-4:end),'.fiff'))
    output = [output '.fif'];
end
if logical(exist(output,'file'))
    curwarn = warning('off', 'backtrace');
    warning('The output file already exist! Now copying...');
    warning(curwarn);
    copyfile(output, [output '.bk']);
end

%-- open & write header
if nargin < 4 || isempty(input)
   fid = fiff_start_file(output);
else
   %--- load input data
   [fid, tree] = fiff_open(input);
   fclose(fid);
    
   timezone         = 5;                        %   Matlab does not know the timezone
   tree.id.secs     = 3600*(24*(now-datenum(1970,1,1,0,0,0))+timezone);
   tree.id.usecs    = 0;                        %   Do not know how we could get this
   
   %--- open output data
   fid = fiff_start_file(output);
   fiff_write_id(fid,FIFF.FIFF_FILE_ID,tree.id)
end

%-- convert unit
if strcmp(unit,'mm')
    convuni = 10^-3;
elseif strcmp(unit,'cm')
    convuni = 10^-2;
else
    convuni = 10^0;
end
trans.trans(1:3,4) = trans.trans(1:3,4) .* convuni;     % apply on move vector

%-- write transform
fiff_write_coord_trans(fid,trans)

%-- write footer
fiff_end_file(fid)

fprintf('Success output tansform-data\n');
