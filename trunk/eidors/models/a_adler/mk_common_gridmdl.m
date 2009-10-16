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
% GREIT V1 (GREIT matrix for reconstruction to circle, 2009)
%   mk_common_gridmdl('GREITc1') - 32x32 with Circle shape
%
% Sheffield Backprojection
%   mk_common_gridmdl('backproj') - 32x32 with Diamond shape
%   mk_common_gridmdl('b2c','backproj') - 32x32 with Circular shape
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
   str= 'b2d';
   RM= 'backproj';
elseif strcmp(str,'GREITc1');
   str= 'b2c';
   RM= 'get_GREIT_c1';
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
if exist('map','var')
   inv_mdl.solve_use_matrix.map = map;
end
inv_mdl.solve = @solve_use_matrix;
inv_mdl.fwd_model.normalize_measurements= 1;
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
   RM = double(GREIT_v10_Circ_Matrix);

function RM = unpack_matrix(packed_matrix);

   % Take a slice
   [x,y]= meshgrid(1:16,1:16);
   ss1 = (y-x)>1 & (y-x)<15;
   sel1 = abs(x-y)>1 & abs(x-y)<15;
   
   [x,y]= meshgrid(-15.5:15.5,-15.5:15.5);
   ss2 = abs(x-y)<25 & abs(x+y)<25 ...
       & x<0 & y<0 & x>=y ;
   sel2 = abs(x-y)<25 & abs(x+y)<25;
 
   % Build up
   BP  = zeros(16^2, 32^2);
   BP(ss1,ss2) = packed_matrix;
   BP  = reshape(BP, 16,16,32,32);

   % Reciprocity
   BP  = BP + permute(BP, [2,1,3,4]);

   % FLIP LR
   el= 16:-1:1;
   BP= BP + BP(el,el,[32:-1:1],:);
   % FLIP UD
   el= [8:-1:1,16:-1:9];
   BP= BP + BP(el,el,:,[32:-1:1]);
   % Transpose
   el= [12:-1:1,16:-1:13];
   BP= BP + permute(BP(el,el,:,:), [1,2,4,3]);

   % Final UD flip to match radiological view (upward toward patient)
   % Here electrodes are connected CW starting from TDC
   BP= BP(:,:,:,[32:-1:1]);

   RM= reshape(BP, 256, [])';
   RM= RM(:,sel1);
% This creates the diamond shape, but we want to leave shape choice later
%  RM= RM(sel2,sel1);


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
    

