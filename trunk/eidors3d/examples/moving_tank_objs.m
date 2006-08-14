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
%   if data_sel is a vector of EIT data, then it will
%      be used for the reconstruction. The first measurement
%      is the homogeneous, and the next are the data
%
% inv_sel
%  2D reconstructions
%   inv_sel = 1   => aa_inv_solve
%   inv_sel = 1.1 => aa_inv_solve (576 elems)
%   inv_sel = 2   => inv_kalman_diff
% 
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.12 2006-08-14 20:48:23 aadler Exp $


if nargin<1; data_sel = 1; end
if nargin<2; inv_sel  = 1; end

if isnumeric(data_sel) & size(data_sel)==[1,1]
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
        load netgen_moving_ball
        vh= vh;
        vi= vi;

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
else
  vh= data_sel(:,1);
  vi= data_sel(:,2:end);
  filename= 'user-data';
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
        imdl.hyperparameter.value= 3e-2;
        time_steps= 1;

        imdl.RtR_prior= @time_smooth_prior;
%       imdl.R_prior= @time_smooth_prior;
        imdl.time_smooth_prior.space_prior= @laplace_image_prior;
        imdl.time_smooth_prior.time_weight= .2;
        imdl.time_smooth_prior.time_steps=  time_steps;

        imdl.solve= @time_prior_solve;
        imdl.time_prior_solve.time_steps=   time_steps;
%       imdl.solve= @np_inv_solve;
%       imdl.solve= 'aa_inv_conj_grad';
end

imgs= inv_solve(imdl,vi,vh);
animate_reconstructions(filename, imgs);
fprintf('file %s.gif created\n',filename);

function rm_rf(dirname)
   if isdir(dirname)==0
       return
   end

   if isunix
       system(['rm -rf "',dirname,'"']);
   else
       system(['rmdir /s /q "',dirname,'"']);
   end
