function moving_tank_objs
% MOVING_TANK_OBJS: create movies of objects moving in tanks
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.3 2005-12-05 16:36:33 aadler Exp $
load guardo_data.mat
%global eidors_colours; eidors_colours.mapped_colour= 255;

imdl= mk_common_model('b2c',16);
imdl.hyperparameter.value= 1e-4;
imdl.RtR_prior.func= 'laplace_image_prior';
imdl= rmfield(imdl,'image_prior');
imdl.solve= 'inv_kalman_diff';

vh= zc_h_demo4;
vi= zc_demo4;
% show_fem(inv_solve( imdl, vi(:,1), vh)); return
imgs= inv_solve(imdl,vi,vh);
mk_movie('back-forth.avi', imgs);
return;
vh= zc_h_demo3;
vi= zc_demo3;

mk_movie('turn.avi', imdl, vi, vh);




% create avi movie fname
% imdl is reconstruction model
% vi is inhomogeneous data sequence
% vh is homogeneous data
function mk_movie(fname, imgs)
   fig=figure;
   set(fig,'DoubleBuffer','on');
   mov = avifile( fname );%, 'Compression', 'RLE' );

   for i=1:length(imgs)
     show_slices(imgs(i));
     F = getframe(gca);
     mov = addframe(mov,F);
   end
   mov = close(mov);
   close(fig);




