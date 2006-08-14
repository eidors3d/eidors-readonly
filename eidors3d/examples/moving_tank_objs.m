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
% $Id: moving_tank_objs.m,v 1.14 2006-08-14 22:43:41 aadler Exp $


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
        vh= homg_tank;
        vi= target_spiral;
        filename= 'netgen-spiral';

    case 31
        load netgen_moving_ball
        vh= homg_tank;
        vi= target_spiral;
        filename= 'netgen-spiral';

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

    case 11
        vv = sub_frame(vi);
        vi = repack_frame(vv);
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';
end

imgs= inv_solve(imdl,vi,vh);
animate_reconstructions(filename, imgs);

function rm_rf(dirname)
   if isdir(dirname)==0
       return
   end

   if isunix
       system(['rm -rf "',dirname,'"']);
   else
       system(['rmdir /s /q "',dirname,'"']);
   end

% simulate condition where frame changes for each
% stimulation. Thus each frame data represents
% a different stimulation pattern
%
% Hardcoded for 16 electrode 2D adjacent stimulation
function vv= sub_frame(vi)
   if ~isstruct(vi);
      viold= vi; vi=struct;
      [st, els]= mk_stim_patterns(16, 1, '{ad}','{ad}', {}, 10);
      for i=1:size(viold,2)
         vi(i).meas= viold(els,i);
      end
   end

   ne= 16; % nelectrodes
   na= ne-3; % data per frame (adjacent)
   for i=1:length(vi)
      f=rem(i-1,ne)+1; %extract jth frame
      vv(i) = eidors_obj('data','', ...
           'configuration', sprintf('Data from %dth stimulation',j), ...
           'meas', vi(i).meas((f-1)*na+(1:na)) );
   end

% native reconstruction will simply push data 
% back together from different frames
%
% Hardcoded for 16 electrode 2D adjacent stimulation

function vv = repack_frame(vi);
   vv= vi;
   l_vi= length(vi);
   ne= 16; % nelectrodes
   na= ne-3; % data per frame (adjacent)
   for i=1:l_vi
      meas= NaN*ones(ne*na,1);
      offset= (1:ne) - min(ne/2,i) + min(l_vi-i-ne/2, 0);
      for j= i+offset;
         f=rem(j-1,ne)+1; %extract jth frame
         meas((f-1)*na+(1:na)) = vi(j).meas;
      end
      vv(i).meas= meas;
   end
