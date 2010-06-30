function VOL = get_elem_volume( fwd_model, map_node )
% GET_ELEM_VOLUME: VOL = get_elem_volume(fwd_model, map_node )
% Calculate volume (or area) of each element in model
%
% If the model has a 'coarse2fine' element, then the
% returned VOL applies to the coarse matrix
%
% if map_node == 1, then calculated volumes are the volume fraction for each node

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin==1; map_node= 0; end

% calculate element volume and surface area
NODE = fwd_model.nodes';
ELEM = fwd_model.elems';
[d,e]= size(ELEM);

VOL=zeros(e,1);
ones_d = ones(1,d);
d1fac = prod( 1:d-1 );
if d > size(NODE,1)
   for i=1:e
       this_elem = NODE(:,ELEM(:,i)); 
       VOL(i)= abs(det([ones_d;this_elem])) / d1fac;
   end
elseif d == 3 % 3D nodes in 2D mesh
   for i=1:e
       this_elem = NODE(:,ELEM(:,i)); 
       d12= det([ones_d;this_elem([1,2],:)])^2;
       d13= det([ones_d;this_elem([1,3],:)])^2;
       d23= det([ones_d;this_elem([2,3],:)])^2;
       VOL(i)= sqrt(d12 + d13 + d23 ) / d1fac;
   end
elseif d == 2 % 3D nodes in 1D mesh (ie resistor mesh)
   for i=1:e
       this_elem = NODE(:,ELEM(:,i)); 
       d12= det([ones_d;this_elem([1],:)])^2;
       d13= det([ones_d;this_elem([2],:)])^2;
       d23= det([ones_d;this_elem([3],:)])^2;
       VOL(i)= sqrt(d12 + d13 + d23 ) / d1fac;
   end
else
   error('mesh size not understood when calculating volumes')
end

if isfield(fwd_model,'coarse2fine')
   VOL= fwd_model.coarse2fine' * VOL;
end

% Calculate the mapping of each element onto the associated node
% Map(i,j) = 1/Ne if elem j has node i
if map_node
   map = sparse( ELEM, ones(d,1)*(1:e), 1/d, size(NODE,2),e);
   VOL = map * VOL;
end
