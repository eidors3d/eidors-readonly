function inv_mdl= mk_common_gridmdl( str, RM)
% MK_COMMON_MODEL: make EIT on reconstruction grids (GREIT)
%
% Usage
%      inv_mdl = mk_common_gridmdl( mdl_string, RecMtx )
%
% 2D models
%   mk_common_gridmdl('b2c', RM)  - 32x32 circular shape, 16 elec
%   mk_common_gridmdl('b2d', RM)  - 32x32 diamond shape, 16 elec
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

if strcmp(str,'backproj')
   str= 'b2d';
   RM= get_Sheffield_Backproj;
elseif strcmp(str,'GREITc1');
   str= 'b2c';
   RM= get_GREIT_c1;
end

n_elec= 16;

switch str
   case 'b2c'
      space = linspace( -1, 1, 32+1 );
      fmdl= mk_grid_model( [], space, space);
      space_avg = conv2(space, [1,1]/2,'valid'); % average of each box
      [x,y] = meshgrid( space_avg, space_avg);
      inside=  (x(:).^2 + y(:).^2)<1.10 ;
      ff = find(~inside);
      fmdl.elems([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine(:,ff)= [];

   case 'b2d'
      space = linspace( -1, 1, 32+1 );
      fmdl= mk_grid_model( [], space, space);
      space_avg = conv2(space, [1,1]/2,'valid'); % average of each box
      [x,y] = meshgrid( space_avg, space_avg);
      inside=  (abs(x(:)) + abs(y(:)))<1.55 ;
      ff = find(~inside);
      fmdl.elems([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine(:,ff)= [];

   otherwise
      error(['mdl_string ',str,' not understood']);
end

fmdl.electrode = mk_electrode_locns( fmdl.nodes, n_elec );

inv_mdl = eidors_obj('inv_model',['mk_common_gridmdl: ',str]);
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
   load GREIT_2009_circ_Matrix.mat
   RM = double(GREITc1);

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
   RM= RM(sel2,sel1);


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
    

