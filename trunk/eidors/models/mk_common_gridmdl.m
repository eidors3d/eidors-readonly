function inv_mdl= mk_common_gridmdl( str, RM)
% MK_COMMON_MODEL: make EIT on reconstruction grids (GREIT)
%
% Usage
%      inv_mdl = mk_common_gridmdl( mdl_string, RecMtx )
%
% 2D models
%   mk_common_gridmdl('b2c', RM)  - 32x32 circular shape, 16 elec
%   mk_common_gridmdl('b2d', RM)  - 32x32 diamond shape, 16 elec
%   mk_common_gridmdl('b2t?', RM) - 32x32 thorax shape, 16 elec
%             - thorax levels 1-5 are provided
%
% Note that the electrodes added to the model are just to 
%   indicate location, it does not necessarily correspond to the
%   Reconstruction Matrix RM provided. 
%
% GREIT MODELS - calculated for library models
%  mk_common_gridmdl('GREIT', PARAMS TO MK_LIBRARY_MODEL)
%    (run mk_library_model('list') to get list)
%  eg. mk_common_gridmdl('GREIT','cylinder_16x1el_coarse')
%  eg. mk_common_gridmdl('GREIT','adult_male_32el')
%
% GREIT V1 (GREIT matrix for reconstruction to circle, 2009)
%   mk_common_gridmdl('GREITc1') - 32x32 with Circle shape
%
% Sheffield Backprojection
%   mk_common_gridmdl('b2d','backproj') - 32x32 with Diamond shape
%   mk_common_gridmdl('b2c','backproj') - 32x32 with Circular shape
%   mk_common_gridmdl('backproj') - 32x32 with Circular shape
%
% COPYRIGHT NOTICE FOR BACKPROJECTION MATRIX:
%   This matrix is copyright DC Barber and BH Brown at
%   University of Sheffield. It may be used free of
%   charge for research and non-commercial purposes.
%   Commercial applications require a licence from the
%   University of Sheffield.
%  

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

name = str;
if strcmp(str,'backproj')
   str= 'b2c';
   RM= 'backproj';
elseif strcmp(str,'GREITc1');
   str= 'b2c';
   RM= 'GREITc1';
else
   name = str;
end

n_elec= 16;

      space = linspace( -1, 1, 32+1 );
      fmdl= mk_grid_model( [], space, space);
      space_avg = conv2(space, [1,1]/2,'valid'); % average of each box
      [x,y] = ndgrid( space_avg, space_avg);
switch str
   case 'b2c'
      inside=  (x(:).^2 + y(:).^2)<1.10 ;

   case 'b2d'
      inside=  (abs(x(:)) + abs(y(:)))<1.55 ;

   case 'b2t1'; inside = inside_thorax(x,y,1);
   case 'b2t2'; inside = inside_thorax(x,y,2);
   case 'b2t3'; inside = inside_thorax(x,y,3);
   case 'b2t4'; inside = inside_thorax(x,y,4);
   case 'b2t5'; inside = inside_thorax(x,y,5);
   
   case 'GREIT'; inv_mdl = greit_library_model(RM); return

   otherwise
      error(['mdl_string ',str,' not understood']);
end

if isstr(RM)
   switch RM
     case 'backproj'; RM= get_Sheffield_Backproj;
     case 'GREITc1';  RM= get_GREIT_c1;
     otherwise;       error('RM string not understood');
   end 
end
      ff = find(~inside);
      fmdl.elems([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine(:,ff)= [];

fmdl.electrode = mk_electrode_locns( fmdl.nodes, n_elec );

inv_mdl = eidors_obj('inv_model',['mk_common_gridmdl: ',name]);
inv_mdl.reconst_type= 'difference';
inv_mdl.fwd_model= fmdl;
inv_mdl.solve_use_matrix.RM = resize_if_reqd(RM,inside);
% solve_use_matrix has a c2f mapping field map, which
% it may be useful to populate in this case
%  inv_mdl.solve_use_matrix.map = map;
inv_mdl.solve = @solve_use_matrix;
inv_mdl.fwd_model = mdl_normalize(inv_mdl.fwd_model, 1);
[st, els]= mk_stim_patterns(16, 1, '{ad}','{ad}', {}, 10);
inv_mdl.fwd_model.meas_select= els;

function inside = inside_thorax(x,y,level);
   [x_bdy, y_bdy ] = thorax_geometry(level,1); % normalized
   inside = inpolygon(x(:), y(:), x_bdy, y_bdy); 

function RM = resize_if_reqd(RM,inside);
   szRM = size(RM,1);
   if sum(inside) == szRM
      % RM is fine
   elseif size(inside,1) == szRM
      RM = RM(inside,:);
   else
      error('mismatch in size of provided RecMatrix');
   end

function RM = get_Sheffield_Backproj
   load Sheffield_Backproj_Matrix.mat
   RM = unpack_matrix(Sheffield_Backproj_Matrix);

function RM= get_GREIT_c1;
   load GREIT_v10_Circ_Matrix.mat
   RM = unpack_matrix(GREIT_v10_Circ_Matrix);


function RM= unpack_matrix(MM)
   RM = unpack_reconst_matrix(MM, 16, 32);


function elec = mk_electrode_locns( nodes, n_elec );
   phi = linspace(0, 2*pi, n_elec+1); 
   phi(end) = [];
   rad = 1;

   for i= 1:n_elec
      posn_x = rad*sin(phi(i));
      posn_y = rad*cos(phi(i));
      dist = (nodes(:,1)-posn_x).^2 + (nodes(:,2)-posn_y).^2;
      [jnk,e_node] = min(dist);

      elec(i).z_contact = 0.001;
      elec(i).nodes     = e_node;
   end
   
function inv_mdl = greit_library_model(str)
fmdl = mk_library_model(str);
fmdl = mdl_normalize(fmdl,1);
nelec = numel(fmdl.electrode);
fmdl.stimulation = ...
   mk_stim_patterns(nelec,1,[0,1],[0,1],{'no_meas_current'}, 1);

opt.noise_figure = 1;
opt.target_size = 0.1;
opt.square_pixels = 1;
if strcmp(str(1:8),'cylinder');
   opt.distr = 0; % seems best, but only works for cylinders
else
   opt.distr = 4; % best bet for the time being
end

inv_mdl = mk_GREIT_model(fmdl,0.25,5,opt);
