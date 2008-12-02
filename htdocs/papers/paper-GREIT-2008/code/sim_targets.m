function sim_targets( savefile )
% simulate radial movement $Id$

   load ng_cyl_mdl.mat
   stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}', {'no_meas_current'}, 1);
   fmdl.stimulation = stim_pat;
   fmdl.elems= double(fmdl.elems);

   imdl= eidors_obj('inv_model','','fwd_model',fmdl);
   imdl.jacobian_bkgnd.value = 1;
   img =  calc_jacobian_bkgnd( imdl );
   vh= fwd_solve( img );
   vh= vh.meas;

   mdl_pts = interp_mesh( fmdl, 2); % 10 per elem
   x= mdl_pts(:,1,:);
   y= mdl_pts(:,2,:);
   z= mdl_pts(:,3,:);

   radius= 0.5*( max(fmdl.nodes(:,1))- min(fmdl.nodes(:,1)) );
   zt    =       max(fmdl.nodes(:,3));
   z0    =       min(fmdl.nodes(:,3));
   rp    = 0.1*radius;
   n_pts = 800;
   c2f    = sparse(size(fmdl.elems,1),n_pts);
   xyzr_pt= zeros(4,n_pts);

   j=1;
   for f_frac = linspace(0,1,n_pts);
%     [xp,yp,zp]= simulation_random(f_frac, radius, z0,zt);
      [xp,yp,zp]= simulation_radmove(f_frac, radius, z0,zt);
      xyzr_pt(:,j) = [xp;yp;zp;rp];
      ff=  (x-xp).^2 + (y-yp).^2 + (z-zp).^2 <= rp^2;
   
      c2f(:,j) =  mean(ff,3); j=j+1;
   end

   fmdl.coarse2fine = c2f;
   J= calc_jacobian( fmdl, img);

   vi = vh*ones(1,n_pts) + J;

   xyzr_pt= diag([-1,1,1,1])*xyzr_pt([2,1,3,4],:)/radius;
   save(savefile, 'vh','vi','xyzr_pt');
   return


   params= [0.9,0.05,0.5,0.5]; %max_posn, targ_rad, z_0, z_t

   [vh,vi,xyzr_pt]= simulate_3d_movement(200, fmdl, params, ...
         @simulation_radmove);
%        @simulation_randctr);
%        @simulation_random);

   %Change: mdl geometry at 90 deg; radius is 15

   save(savefile, 'vh','vi','xyzr_pt');

function vh= homog_sim( fmdl );

function [xp,yp,zp]= simulation_randctr(f_frac, radius, z0,zt);
   rad= rand(1)*radius;
   phi= 2*pi*rand(1);
   xp = rad*sin(phi);
   yp = rad*cos(phi);
%  zp = mean([zt,z0]) + 0.05*std([zt,z0])*randn(1);
   zp = mean([zt,z0]);


function [xp,yp,zp]= simulation_random(f_frac, radius, z0,zt);
   lim = (radius*(1-0.05-0.05))^2;
   while 1
     xp = radius*2*(rand(1)-0.5);
     yp = radius*2*(rand(1)-0.5);
     if xp^2 + yp^2 < lim; break;end
   end
   zp = mean([zt,z0]) + 0.05*std([zt,z0])*randn(1);

function [xp,yp,zp]= simulation_radmove(f_frac, radius, z0,zt);
   rp= f_frac*radius; 
   cv= 2*pi*f_frac * 73;
   xp= rp * cos(cv);
   yp= rp * sin(cv);
   zp = mean([zt,z0]);

