function moving_tank_objs
% MOVING_TANK_OBJS: create movies of objects moving in tanks
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.1 2005-12-05 11:54:20 aadler Exp $
load guardo_data.mat

imdl= mk_common_model('c2c',16);
imdl.hyperparameter.value= 1e-4;

vh= zc_h_demo3;
vi= zc_demo3;

mk_movie('turn.avi', imdl, vi, vh);



% create avi movie fname
% imdl is reconstruction model
% vi is inhomogeneous data sequence
% vh is homogeneous data
function mk_movie(fname, imdl, vi, vh)
   fl= size(vi,2);

   fig=figure;
   set(fig,'DoubleBuffer','on');
   mov = avifile( fname );

   for i=1:fl
     show_slices(inv_solve(imdl,vi(:,i),vh));
     F = getframe(gca);
     mov = addframe(mov,F);
   end
   mov = close(mov);
   close(fig);




