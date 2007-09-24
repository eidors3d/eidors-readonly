% Map elements onto mesh
% $Id: moving_ball02.m,v 1.3 2007-09-24 17:56:59 aadler Exp $

% Load models at start

load ng_mdl_16x1_ptrs; % calculated from mk_mesh_sample_array
xx=xyz(:,1); yy=xyz(:,2); zz=xyz(:,3);

stimulation= mk_stim_patterns(16,1,'{ad}','{ad}', ...
                      {'meas_current', 'do_redundant'});

for mdlidx= 1:3
   if     mdlidx==1; fmdl= ng_mdl_16x1_coarse;
   elseif mdlidx==2; fmdl= ng_mdl_16x1_fine;
   elseif mdlidx==3; fmdl= ng_mdl_16x1_vfine;
   end

   fmdl.stimulation= stimulation;
 
   f_no = 200;              % number to simulate
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
       xyzr_pt(:,i)= [xp;-yp;zp;rp]; % -y because images and axes are reversed

       ff= find( (xx-xp).^2 + (yy-yp).^2 + (zz-zp).^2 <= rp^2 )';
       obj_n= sparse( eptr(ff,mdlidx),1,1, n_elems, 1);

       img.elem_data= 1 + target_conductivity * (obj_n./eptr_n);

       vi(i)= fwd_solve( img );
%show_fem(img);drawnow; keyboard
   end

   % convert to data matrix
   vi= [vi(:).meas]; 
   vh= [vh(:).meas];
   if     mdlidx==1; save move_ball_16x1_adj_coarse.mat vi vh xyzr_pt
   elseif mdlidx==2; save move_ball_16x1_adj_fine.mat   vi vh xyzr_pt
   elseif mdlidx==3; save move_ball_16x1_adj_vfine.mat  vi vh xyzr_pt
   end
   clear vi vh;
end
