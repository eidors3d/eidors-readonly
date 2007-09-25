% Map elements onto mesh
% $Id: moving_ball_16x2_02.m,v 1.1 2007-09-25 11:21:44 aadler Exp $

% Load models
fmdl= ng_mdl_16x2_vfine;
fmdl= ng_mdl_16x2_fine;
fmdl= ng_mdl_16x2_coarse;

mdlidx= 1;
load ng_mdl_16x2_ptrs; % calculated from mk_mesh_sample_array
xx=xyz(:,1); yy=xyz(:,2); zz=xyz(:,3);
 
   f_no = 2  ;              % number to simulate
   target_conductivity= .2; % Delta coneductivity for contrast

   xyzr_pt= [];

   tank_maxdims = max(fmdl.nodes) - min(fmdl.nodes);
   tank_radius = 0.5*min(tank_maxdims(1:2));
   tank_height = tank_maxdims(3);

   rp= .05*tank_radius; 
   path_radius= 2/3*tank_radius;

   % Homogeneous model
   n_elems= size(fmdl.elems,1);
   img= eidors_obj('image','simulate_movement', ...
                   'fwd_model', fmdl, ...
                   'elem_data', ones(n_elems,1) );
   vh= fwd_solve(img);

   % n pointers in each element
   ff= ~isnan(eptr(:,mdlidx));
   eptr_n= sparse( eptr(ff,mdlidx),1,1, n_elems, 1);

   for i=1:f_no
       f_frac= (i-1)/f_no;
       fprintf('simulating %d / %d \n',i,f_no);

       xp= path_radius * cos(f_frac*2*pi);
       yp= path_radius * sin(f_frac*2*pi);
       zp= tank_height / 2;
       xyrz_pt(:,i)= [xp;-yp;zp;rp]; % -y because images and axes are reversed

       ff= find( (xx-xp).^2 + (yy-yp).^2 + (zz-zp).^2 <= rp^2 )';
       obj_n= sparse( eptr(ff,mdlidx),1,1, n_elems, 1);

       img.elem_data= 1 + target_conductivity * (obj_n./eptr_n);

       vi(i)= fwd_solve( img );
%show_fem(img);drawnow; keyboard
    end

% convert to data matrix
vi= [vi(:).meas]; 
vh= [vh(:).meas];
