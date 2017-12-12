function VOL = get_elem_volume( fwd_model, map_node )
% GET_ELEM_VOLUME: VOL = get_elem_volume(fwd_model, map_node )
% Calculate volume (or area) of each element in model
%
% If the model has a 'coarse2fine' element, then the
% returned VOL applies to the coarse matrix (unless map_node <0)
%
% if map_node < 0, do not apply coarse2fine (if it exists)
%
% if abs(map_node) == 1, then calculated volumes are the volume fraction for each node
% BUGS: can't currently apply coarse2fine on map_node.

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

% If it looks like a fwd model
if ~isfield(fwd_model,'type');
  if isfield(fwd_model,'nodes') && isfield(fwd_model,'elems')
     fwd_model.type = 'fwd_model';
  end
end
   

switch fwd_model.type
  case 'fwd_model'; % do nothing, we're ok
  case 'rec_model'; % do nothing, we're ok
  case 'inv_model'; fwd_model = fwd_model.fwd_model;
  case 'image';     fwd_model = fwd_model.fwd_model;
  otherwise;
     error('get_elem_volume: expecting fwd_model, got %s', ...
           fwd_model.type)
end

if nargin==1; map_node= 0; end

% calculate element volume and surface area
NODE = fwd_model.nodes';
ELEM = fwd_model.elems';

copt.cache_obj = {NODE, ELEM,map_node};
copt.fstr = 'elem_volume';
copt.boost_priority = -4;

VOL = eidors_cache(@calc_volume, {NODE, ELEM}, copt);
if isfield(fwd_model,'coarse2fine') && map_node>=0
   VOL= fwd_model.coarse2fine' * VOL;
end

% Calculate the mapping of each element onto the associated node
% Map(i,j) = 1/Ne if elem j has node i
if abs(map_node)==1
   [d,e]= size(ELEM);
   map = sparse( ELEM, ones(d,1)*(1:e), 1/d, size(NODE,2),e);
   VOL = map * VOL;
end

function VOL = calc_volume(NODE,ELEM)
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
    


function do_unit_test
  imdl = mk_common_model('a2c2',8);
  tt = 0.03125*ones(4,1);
  out = get_elem_volume(imdl.fwd_model);
  unit_test_cmp('fmdl:',  out(1:4), tt, 1e-10);
  out = get_elem_volume(imdl);
  unit_test_cmp('imdl:',  out(1:4), tt, 1e-10);
  out = get_elem_volume(mk_image(imdl,1));
  unit_test_cmp('image:', out(1:4), tt, 1e-10);

  fmdl = imdl.fwd_model;
  unit_test_cmp('fmdl (size 1):',  size(out,1), [num_elems(fmdl)]);
  out0= get_elem_volume(fmdl,0);
  unit_test_cmp('fmdl (size 1):',  size(out0,1), [num_elems(fmdl)]);

  out1= get_elem_volume(fmdl,1);
  unit_test_cmp('fmdl (size 1):',  size(out1,1),[num_nodes(fmdl)]);
  unit_test_cmp('fmdl (nodes 1):',  sum(out), sum(out1), 1e-10);
  ttn(1) = 0.041666666666667; ttn(2:5) = 0.088388347648318;
  unit_test_cmp('fmdl (nodes 2):',  out1(1:5), ttn(:), 1e-10);

  fmdl.coarse2fine = zeros(num_elems(fmdl),2);
  fmdl.coarse2fine(1:4,1) = 1;
  fmdl.coarse2fine(5:end,2) = 1;
  outc= get_elem_volume(fmdl,-2); % don't use
  unit_test_cmp('fmdl (c2f 1):',  out0, outc);
  outc= get_elem_volume(fmdl,0); % use c2f
  unit_test_cmp('fmdl (c2f 2):',  outc(1), sum(tt), 1e-10);
  unit_test_cmp('fmdl (c2f 3):',  sum(outc), sum(out0), 1e-10);
  unit_test_cmp('fmdl (c2f 4):',  size(outc), [2,1]);

  outc= get_elem_volume(fmdl,-1); % use c2f and map_node
  unit_test_cmp('fmdl (nodes 2):',  outc(1:5), ttn(:), 1e-10);
