function val = img_point_mapper(img, pts, maptype )
%IMG_POINT_MAPPER - image values at points
% V = IMG_POINT_MAPPER(IMG, P, MAPTYPE ) returns the matrix of values at
% points P (N x 3). 
%   IMG     - EIDORS image structure
%   P       - a list of points (N x 3)
%   MAPTYPE - specifies the value returned for each point:
%     'elem'      - value of the enclosing elements
%     'node'      - value of the nearest node
%     'nodeinterp'- interpoloation of node values of the enclosing element
%   IMG_POINT_MAPPER uses GET_IMG_DATA to obtain image values. Data must be
%   defined on elements or nodes as appropriate for the chosen MAPTYPE.
%   If MAPTYPE is not specified, 'elem' or 'nodeinterp' is selected.
%
% NOTE: Only works in Matlab.
%
% See also: GET_IMG_DATA, MDL_SLICE_MAPPER

% TODO: 
%   - support for coarse2fine

% (C) 2021 Bartek Grychtol. License: GPL version 2 or version 3
% $Id$

data = get_img_data(img);
fmdl = img.fwd_model;
if nargin < 3
    if size(data,1) == size(img.fwd_model.nodes,1)
        maptype = 'nodeinterp';
    else
        maptype = 'elem';
    end
end

TR = eidors_cache(@triangulation,{fmdl.elems, fmdl.nodes});

switch maptype
    case 'elem'
        id = pointLocation(TR, pts);
        val = data(id,:);
    case 'node'
        id = nearestNeighbor(TR, pts);
        val = data(id,:);
    case 'nodeinterp'
        [id, bc] = pointLocation(TR, pts);
        n_pts = size(pts,1);
        map = builtin('sparse', repelem((1:n_pts)',1,4), fmdl.elems(id,:), bc, n_pts, size(fmdl.nodes,1));
        val = map * data;
    otherwise
        error('Maptype must be ''elem'', ''node'', or ''nodeinterp''.')
end
