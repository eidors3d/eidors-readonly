function imgr= moving_tank_objs(data_sel, inv_sel)
% MOVING_TANK_OBJS: create movies of objects moving in tanks
% Usage:
% imgr= moving_tank_objs(data_sel, invmdl)
%
% imgr = reconstructed images
%
% data_sel select data_sel to use
%   data_sel: 10 => target across tank data
%             20 => target around tank data
%             21 => target around tank data (10 positions)
%             30 => from netgen simulations
%             40 => show netgen FEM simulation
%
% inv_sel
%  2D reconstructions
%   inv_sel = 1   => aa_inv_solve
%   inv_sel = 1.1 => aa_inv_solve (576 elems)
%   inv_sel = 2   => inv_kalman_diff
% 
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.11 2006-08-12 00:01:24 aadler Exp $


if nargin<1; data_sel = 1; end
if nargin<2; inv_sel  = 1; end

switch data_sel
    case 10
        load montreal_data_1995;
        vh= zc_h_demo4;
        vi= zc_demo4;
        filename= 'target-across';

    case 20
        load montreal_data_1995;
        vh= zc_h_demo3;
        vi= zc_demo3;
        filename= 'target-around';

    case 21
        load montreal_data_1995;
        vh= zc_h_demo3;
        vi= zc_demo3(:,1:2:20);
        filename= 'target-around10';

    case 30
        load all_models_for_moving_ball
        stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}',{'no_meas_current'},1);
        vh= do_simulation( tank_img_homg, stim_pat);
        for i=1:2 %length(tank_img)
            vi(:,i)= do_simulation( tank_img(i), stim_pat);
        end
        filename= 'ball-model';
        % FIXME: something weird is happening here

    case 40
        dirname= 'tmp_mk_movie2';
        rm_rf( dirname );
        mkdir( dirname );

        load tank_mdls
        for i=1:length(tank_mdls);
            show_fem( tank_mdls(i) );
            print('-dpng', ...
                sprintf('%s/mdl%05d.png',dirname, i) );
        end
        fname= 'move_ball';

        ld_lib_path= sys_dep;
        retval= system(sprintf( ...
            '%s convert -delay 25 %s/img*.png -loop 0 %s.gif', ...
            ld_lib_path, dirname, fname ));
        if retval~=0
            error('please ensure the imagemagick convert program is in your path');
        end
        rm_rf(dirname);

    otherwise
        error('Don''t recognize data_sel');
end

switch inv_sel
    case 1
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';

    case 1.1
        imdl= mk_common_model('c2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';

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

    case 4
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        time_steps= 1;

        imdl.RtR_prior= @time_smooth_prior;
        imdl.time_smooth_prior.space_prior= @laplace_image_prior;
        imdl.time_smooth_prior.time_weight= .5;
        imdl.time_smooth_prior.time_steps=  time_steps;

        imdl.solve= @time_prior_solve;
        imdl.time_prior_solve.time_steps=   time_steps;
%       imdl.solve= @np_inv_solve;
%       imdl.solve= 'aa_inv_conj_grad';
end

imgs= inv_solve(imdl,vi,vh);
mk_movie2(filename, imgs);
fprintf('file %s.gif created\n',filename);



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

   r_img= show_slices(imgs);
   c_img = calc_colours( r_img);
   out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
   cmap= colormap;

   for i=1:length(imgs)
     imwrite(out_img(:,:,i),cmap, ...
            sprintf('%s/img%05d.png',dirname, i), 'png');
   end

   ld_lib_path= sys_dep;

   retval= system(sprintf( ...
       '%s convert -delay 25 %s/img*.png -loop 0 %s.gif', ...
       ld_lib_path, dirname, fname ));
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

% work around stupid matlab bugs
function ld_lib_path= sys_dep;
   ld_lib_path='';
   if  strfind(system_dependent('getos'),'Linux')
     s=ver;
      if str2num(s.Version)>=7
        %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks external progs
          ld_lib_path='LD_LIBRARY_PATH=;';
      end      
   end    
