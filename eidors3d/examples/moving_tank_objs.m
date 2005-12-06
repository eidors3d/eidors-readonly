function moving_tank_objs
% MOVING_TANK_OBJS: create movies of objects moving in tanks
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.4 2005-12-06 00:16:41 aadler Exp $
load guardo_data.mat
%global eidors_colours; eidors_colours.mapped_colour= 255;

imdl= mk_common_model('b2c',16);
imdl.hyperparameter.value= 1e-2;
imdl.RtR_prior.func= 'laplace_image_prior';
if 0
    imdl.solve= 'inv_kalman_diff';
else
    imdl.solve= 'aa_inv_solve';
end

vh= zc_h_demo4;
vi= zc_demo4;
imgs= inv_solve(imdl,vi,vh);
%mk_movie('back-forth.gif', imgs);
mk_movie2('back-forth.gif', imgs);

vh= zc_h_demo3;
vi= zc_demo3;
imgs= inv_solve(imdl,vi,vh);

%mk_movie('turn.avi', imgs);
mk_movie2('turn.gif', imgs);




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

function mk_movie2(fname, imgs)
   mkdir('tmp_mk_movie2');
   cmap= colormap;
   for i=1:length(imgs)
     img= show_slices(imgs(i));
     imwrite(img,cmap, ...
            sprintf('tmp_mk_movie2/img%05d.png',i), 'png');
   end
   retval= system(sprintf( ...
       'convert -delay 25 tmp_mk_movie2/img*.png -loop 0 %s', ...
       fname ));
   if retval~=0
       error('please ensure the imagemagick convert program is in your path');
   end
   if isunix
       system('rm -r tmp_mk_movie2');
   else
       system('rm -r -d tmp_mk_movie2');
   end



