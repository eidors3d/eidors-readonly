function imgr= moving_tank_objs(data_sel, inv_sel)
% MOVING_TANK_OBJS: create movies of objects moving in tanks
% Usage:
% imgr= moving_tank_objs(data_sel, invmdl)
%
% imgr = reconstructed images
%
% data_sel select data_sel to use
%   data_sel = 1 => target across tank data
%   data_sel = 2 => target around tank data
%
% inv_sel
%   inv_sel = 1 => 2D reconstruction=> aa_inv_solve
%   inv_sel = 2 => 2D reconstruction=> inv_kalman_diff
% 
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.6 2006-01-24 02:53:26 aadler Exp $


if nargin<1; data_sel = 1; end
if nargin<2; inv_sel  = 1; end

switch data_sel
    case 1
        load montreal_data_1995;
        vh= zc_h_demo4;
        vi= zc_demo4;
        filename= 'target-across';

    case 2
        load montreal_data_1995;
        vh= zc_h_demo3;
        vi= zc_demo3;
        filename= 'target-around';

    case 3
        load all_models_for_moving_ball
        stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}',{'no_meas_current'},1);
        vh= do_simulation( tank_img_homg, stim_pat);
        for i=1:2 %length(tank_img)
            vi(:,i)= do_simulation( tank_img(i), stim_pat);
        end
        filename= 'ball-model';
        % FIXME: something weird is happening here

    otherwise
        error('Don''t recognize data_sel');
end

switch inv_sel
    case 1
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'aa_inv_solve';

    case 1.1
        imdl= mk_common_model('c2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'aa_inv_solve';

    case 2
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'inv_kalman_diff';

    case 3
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'aa_inv_conj_grad';
end

imgs= inv_solve(imdl,vi,vh);
mk_movie2(filename, imgs);



% create avi movie fname
% imdl is reconstruction model
% vi is inhomogeneous data sequence
% vh is homogeneous data
function mk_movie(fname, imgs)
   fig=figure;
   set(fig,'DoubleBuffer','on');
   mov = avifile( [fname ,'.avi'] );%, 'Compression', 'RLE' );

   for i=1:length(imgs)
     show_slices(imgs(i));
     F = getframe(gca);
     mov = addframe(mov,F);
   end
   mov = close(mov);
   close(fig);

function mk_movie2(fname, imgs)
   calc_colours('mapped_colour', 127);
   dirname= 'tmp_mk_movie2';
   rm_rf( dirname );
   mkdir( dirname );
   for i=1:length(imgs)
     img= show_slices(imgs(i));
     cmap= colormap;
     imwrite(img,cmap, ...
            sprintf('%s/img%05d.png',dirname, i), 'png');
   end
   retval= system(sprintf( ...
       'convert -delay 25 %s/img*.png -loop 0 %s.gif', ...
       dirname, fname ));
   if retval~=0
       error('please ensure the imagemagick convert program is in your path');
   end
   rm_rf(dirname);

function rm_rf(dirname)
   if isdir(dirname)==0
       return
   end

   if isunix
       system(['rm -rf "',dirname,'"']);
   else
       system(['rmdir /s /q "',dirname,'"']);
   end


% remove all electrodes except the first 16.
% apply stimulation pattern and return data
function vv= do_simulation( img, stim_pat)

   img.fwd_model.electrode= ...
      img.fwd_model.electrode(1:16);

   img.fwd_model.stimulation= stim_pat;

   vv= fwd_solve( img);
