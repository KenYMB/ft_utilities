function [Center, Radius] = ft_neurofif_pick_dig_for_fitsphere(input, output, ignorepts)
%
% [center, radius] = FT_NEUROFIF_PICK_DIG_FOR_FITSPHERE(input, output [,ignore_points])
%
% Load isotrak points from 'input', and pick-up appropriate points to
% estimate sphere model.
% If second input argument exists, the appropriate points are saved as 'output'.
% If 'ignore_points' is input, the indicated index of extra points are
% considered as inappropriate points.
%
% -- output --
%   center  = 1x3-array,    the center of the fitted sphere.  [mm]
%   radius  = number,       the radius of the fitted sphere.  [mm]
% 
% -- input --
%   input   = strings,      filename to load isotrak points.
%   output  = strings,      filename to write out picked-up isotrak points.
%                           The positions are saved in head coordinate.  [m]
%                           If the filename extension is '.fif' or '.fiff',
%                           the data is output as neuromag-fif file.
%                           Otherwise, the positions of extra points are
%                           saved as text file.
%   ignorepts = 1xN-array,  the index of extra points to remove.
% 
% -- example --
% center = ft_neurofif_pick_dig_for_fitsphere('MEG_raw.fif');
% fid = fopen('center_sphere','w');
% fprintf(fid,'%f %f %f', center);                  % save as 'mm'
% fclose(fid);
% 
% -- example for UNIX --
% ft_neurofif_pick_dig_for_fitsphere('MEG_raw.fif','headpoints');
% cmd_fit   = sprintf('/neuro/bin/util/fit_sphere_to_points %s','headpoints');
% [status, spherefit] = unix(cmd_fit);
% spherefit = str2num(spherefit)*10^3;              % convert into 'mm'
% center    = spherefit(1:3);
% radius    = spherefit(4);
% 

% using: fieldtrip(, ft_hastoolbox, mne), ft_private(, fitsphere)

% 20170815 Yuasa

me='original:ft_neurofif_pick_dig_for_fitsphere';

narginchk(1,3)

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

isfif  = false;
if nargin < 2 || isempty(output)
    output = '';
elseif length(output)>=5 && (strcmp(output(end-3:end),'.fif') || strcmp(output(end-4:end),'.fiff'))
    isfif = true;
end
if nargin < 3
    ignorepts = [];
end
if logical(exist(output,'file'))
    curwarn = warning('off', 'backtrace');
    warning('The output file already exist! Now copying...');
    warning(curwarn);
    copyfile(output, [output '.bk']);
end
%
%   Open the file
%
fprintf(1,'Reading %s ...',input);
[ fid, tree ] = fiff_open(input);
%
%   Locate the data of interest
%
diginfo = fiff_dir_tree_find(tree,FIFF.FIFFB_ISOTRAK);
if isempty(diginfo)
    fclose(fid);
    error(me,'Could not find insotrak data');
end
%
%   Load the data of interest
%
clear digpoints
for ilp = 1:length(diginfo.dir)
    digpoints(ilp) = fiff_read_tag(fid,diginfo.dir(ilp).pos);
end
fclose(fid);
%
%
%   Identify the digitizer position
%
digheader = digpoints([digpoints.kind]~=FIFF.FIFF_DIG_POINT);     % bakup other info
digpoints([digpoints.kind]~=FIFF.FIFF_DIG_POINT) = [];
digpoints = [digpoints.data];
extralist = [digpoints.kind]==FIFF.FIFFV_POINT_EXTRA;
extra     = [digpoints(extralist).r];
%
%   Detect invalid points
%
validlist = true(1,size(extra,2));
validlist(ignorepts) = false;                                   %-- check ignore points

validlist = validlist & extra(3,:)>=0;                          %-- check if z>=0

distlist  = zeros(size(validlist));                             %-- check the minimum distance to detect miss points
for ilp=1:length(distlist)
    tmpdist = rms(repmat(extra(:,ilp),1,length(distlist)) - extra,1);
    tmpdist(ilp) = [];  tmpdist(find(tmpdist==min(tmpdist),1)) = [];
    distlist(ilp) = min(tmpdist);       % 2nd min
end
validlist = validlist & (distlist < 0.02);                      %   check 2nd min distance < 2cm
%
%   Remove invalid points
%
extra = double(extra(:,validlist)');                            % trans & convert into double
extralist(extralist) = ~validlist;                              % invalid list
digpoints(extralist) = [];
extralist  = [digpoints.kind]==FIFF.FIFFV_POINT_EXTRA;
extraident = num2cell(1:length(find(extralist)));
[digpoints(extralist).ident] = extraident{:};
%
%   Output variables
%
try
    ft_private('silent');
      [Center, Radius] = fitsphere(extra);
    ft_private('silent','remove');
    Center = Center .* 10^3;                                    % convert m -> mm
    Radius = Radius .* 10^3;                                    % convert m -> mm
catch
    fprintf(2,'\nfail to find the function to fit sphere!\n');
    Center = [];    Radius = [];
end
%
%   Save output file
%
if ~isempty(output)
  fprintf(1,'Writing %s ...',output);
  if isfif          % save as fif
    ft_write_dig(output, digheader, digpoints);
  else              % save as txt
    save(output, 'extra', '-ascii');
  end
  fprintf('\nSuccess output isotrak-data\n');
else
  fprintf('Finish\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set to ouput (instead of 'fiff_write_dig_file')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ft_write_dig(output, digheader, digpoints)

global FIFF;
%
%   Start writing...
%
fid  = fiff_start_file(output);
fiff_start_block(fid,FIFF.FIFFB_ISOTRAK);
%
% Ready to start
%
for ilp = 1:length(digheader)
    fiff_write_int(fid,digheader(ilp).kind,digheader(ilp).data);
end
%
% Write in each
%
for ilp = 1:length(digpoints)
    fiff_write_dig_point(fid,digpoints(ilp));
end
%
%   Finish writing...
%
fiff_end_block(fid,FIFF.FIFFB_ISOTRAK);
fiff_end_file(fid);

