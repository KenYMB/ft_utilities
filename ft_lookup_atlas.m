function [varargout] = ft_lookup_atlas(atlas, pos, varargin)

% ATLAS_LOOKUP determines the anatomical label of a location in the given atlas.
%
% Use as
%   label = atlas_lookup(atlas, pos, ...);
% 
%   labels are sorted along with the distance 
%   inputcoords between atlas and pos have to be same
% 
% If atlas includes several labelsets, 
%   [label1, label2, ...] = atlas_lookup(atlas, pos, ...);
% 
%   outputs label names for each labelsets.
%
% Optinal input arguments should come in key-value pairs and can include
%   'queryrange'   = number, should be 1, 3, 5, 7, 9, 11 or 13 (default = 5)
%                    (1=just the position, 5=+-2 range from the position)
%   'outputs'      = maximum number of output atlas (default = [])
% 
% Dependent on the input coordinates and the coordinates of the atlas, the
% input positions are transformed betweem MNI and Talairach-Tournoux coordinates.
% See http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml for more details.

% Copyright (C) 2005-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% 20170914 Yuasa: based on atlas_lookup

% get the optional input arguments
queryrange  = ft_getopt(varargin, 'queryrange', 5);
outputs     = ft_getopt(varargin, 'outputs');

if isempty(intersect(queryrange, [1 3 5 7 9 11 13]))
  error('incorrect query range, should be one of [1 3 5 7 9 11 13]');
end

if size(pos,1)==3 && size(pos,2)~=3
  % transpose the input positions to get Nx3
  pos = pos';
end

% determine which field(s) to use to look up the labels,
% and whether these are boolean or indexed
fn = fieldnames(atlas);
isboolean = false(numel(fn),1);
isindexed = false(numel(fn),1);
for i=1:length(fn)
  if islogical(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
    isboolean(i) = true;
  elseif isnumeric(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
    isindexed(i) = true;
  end
end
if any(isindexed)
  % let the indexed prevail
  fn = fn(isindexed);
  isindexed = 1;
elseif any(isboolean)
  % use the boolean
  fn = fn(isboolean);
  isindexed = 0;
end

num = size(pos,1);
sel = cell(1,numel(fn));
dis = cell(1,numel(fn));
varargout = cell(1,numel(fn));

outputnum = zeros(numel(fn),1);
for k=1:numel(fn)
    labelnum     = length(atlas.([fn{k} 'label']));
    if isempty(outputs),        outputnum(k) = labelnum;
    elseif length(outputs)<k,   outputnum(k) = labelnum;
    else                        outputnum(k) = min(outputs(k),labelnum);
    end
    varargout{k} = repmat({''},num,outputnum(k));
end

% convert the atlas head coordinates into voxel coordinates
vox  = ft_warp_apply(inv(atlas.transform), pos);

for i=1:num

  % this is the center voxel
  ijk_center = vox(i,:);

  if isindexed

    [di,dj,dk] = ndgrid((-(queryrange-1)/2):1:((queryrange-1)/2), ...
                        (-(queryrange-1)/2):1:((queryrange-1)/2), ...
                        (-(queryrange-1)/2):1:((queryrange-1)/2));
    dijk = [di(:),dj(:),dk(:)];
    
    % search in a cube around the center voxel
    ijk = repmat(round(ijk_center),queryrange.^3,1) + dijk;
    ijk(any([ijk<1, ijk(:,1)>atlas.dim(1), ...
                    ijk(:,2)>atlas.dim(2), ...
                    ijk(:,3)>atlas.dim(3)], 2), :) = [];
    ijk_ind  = sub2ind(atlas.dim,ijk(:,1),ijk(:,2),ijk(:,3));

    for k=1:numel(fn)
      sel{k} = atlas.(fn{k})(ijk_ind);
      dis{k} = sum((ijk - repmat(ijk_center,length(ijk_ind),1)).^2,2);
    end
    
    for k = 1:numel(fn)
      if ~isempty(sel{k})
        % sort along with distance
        seltmp = sortrows(cat(2,dis{k},reshape(double(sel{k}),[],1)));
        labidx = uint16(setdiff(unique(seltmp(:,2),'stable'),0,'stable'));
        if length(labidx) > outputnum(k)
            labidx = labidx(1:outputnum(k));
        end
        varargout{k}(i,1:length(labidx)) = atlas.([fn{k} 'label'])(labidx);
      end
    end
  else
    error('support for atlases that have a probabilistic segmentationstyle is not supported yet');
  end
end

%label = unique(atlas.descr.name(sel));

