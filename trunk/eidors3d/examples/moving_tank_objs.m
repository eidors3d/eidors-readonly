function imgs= moving_tank_objs(data_sel, inv_sel, options)
% MOVING_TANK_OBJS: create movies of objects moving in tanks
% Usage:
% imgs= moving_tank_objs(data_sel, inv_sel, options)
%
% imgs = reconstructed images
%
% data_sel select data_sel to use
%   data_sel: 10 => target across tank data
%             20 => target around tank data
%             21 => target around tank data (10 positions)
%             30 => from netgen simulations (spiral)
%             32 => from netgen simulations (spirograph)
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
% options (not changed if value is NaN);
%   options(1) - hyperparameter
%   options(2) - time_steps
%   options(3) - time_weight
% 
% Create moving objects and tanks
% $Id: moving_tank_objs.m,v 1.16 2006-08-25 00:14:37 aadler Exp $

clim= [];

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
        vi= target_spiral(40:60);
        filename= 'netgen-spiral-part';

    case 32
        load netgen_moving_ball
        vh= homg_tank;
        vi= target_spirograph_slow;
        filename= 'netgen-spirograph-slow';

    case 33
        load netgen_moving_ball
        vh= homg_tank;
        vi= target_spirograph_slow([31:52,72:90]);
        randn('seed',20);
        for i=1:length(vi)
           vi(i).meas = vi(i).meas + 25e-6*randn(size(vi(i).meas));
        end
        filename= 'netgen-spirograph-slow-part';

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

    case 1.2
        imdl= mk_common_model('d2c',16);
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
        time_steps= 0;

        imdl.RtR_prior= @time_smooth_prior;
        imdl.time_smooth_prior.space_prior= @laplace_image_prior;
        imdl.time_smooth_prior.time_weight= .04;
        imdl.time_smooth_prior.time_steps=  time_steps;

        imdl.solve= @time_prior_solve;
        imdl.time_prior_solve.time_steps=   time_steps;

    case 4.1
        imdl= mk_common_model('c2c',16);
        imdl.hyperparameter.value= 3e-1;
        time_steps= 3;

        imdl.RtR_prior= @time_smooth_prior;
        imdl.time_smooth_prior.space_prior= @noser_image_prior;
        imdl.noser_image_prior.exponent= .5;
        imdl.time_smooth_prior.time_weight= .5;
        imdl.time_smooth_prior.time_steps=  time_steps;

        imdl.solve= @time_prior_solve;
        imdl.time_prior_solve.time_steps=   time_steps;

    case 4.2
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 3e-1;
        time_steps= 0;

        imdl.RtR_prior= @time_smooth_prior;
        imdl.time_smooth_prior.space_prior= @tikhonov_image_prior;
        imdl.time_smooth_prior.time_weight= .01;
        imdl.time_smooth_prior.time_steps=  time_steps;

        imdl.solve= @time_prior_solve;
        imdl.time_prior_solve.time_steps=   time_steps;

    case 11
        vv = sub_frame(vi);
        vi = repack_frame(vv);
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';

    case 11.1
        vv = sub_frame(vi);
        vi = repack_frame(vv);
        imdl= mk_common_model('c2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';

    case 12
        vi = sub_frame(vi,vh);
        vh = zeros(length(vi(1).meas),1);
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= @laplace_image_prior;
        imdl.solve= @inv_kalman_diff;
        for i=1:16;
           imdl.inv_kalman_diff.sequence(i).meas_no= (i-1)*13+(1:13);
        end
        imdl.fwd_model = rmfield(imdl.fwd_model,'meas_select');
        imdl.fwd_model.normalize=0;

    otherwise
        error(['inv_sel (' inv_sel ') not recognized']);
end

% options
if nargin>=3
   if length(options) >=1 if ~isnan(options(1))
         imdl.hyperparameter.value= options(1);
         filename = sprintf('%s-hp=%g', filename, options(1));
   end; end
   if length(options) >=2 if ~isnan(options(2))
      imdl.time_smooth_prior.time_steps=  options(2);
      imdl.time_prior_solve.time_steps=   options(2);
      filename = sprintf('%s-ts=%d', filename, options(2));
   end; end
   if length(options) >=3 if ~isnan(options(3))
      imdl.time_smooth_prior.time_weight=  options(3);
      filename = sprintf('%s-tw=%4.2f', filename, options(3));
   end; end
end

imgs= inv_solve(imdl,vi,vh);
animate_reconstructions(filename, imgs, clim);

% simulate condition where frame changes for each
% stimulation. Thus each frame data represents
% a different stimulation pattern
%
% if both vi and vh are provided, calculate difference
%  signal
%
% Hardcoded for 16 electrode 2D adjacent stimulation
function vv= sub_frame(vi,vh)
   if ~isstruct(vi);
      viold= vi; vi=struct;
      [st, els]= mk_stim_patterns(16, 1, '{ad}','{ad}', {}, 10);
      for i=1:size(viold,2)
         vi(i).meas= viold(els,i);
      end
      if nargin==2
         vh = struct('meas', vh(els,1));
      end
   end

   ne= 16; % nelectrodes
   na= ne-3; % data per frame (adjacent)
   for i=1:length(vi)
      f=rem(i-1,ne)+1; %extract jth frame
      idx= (f-1)*na+( 1:na );
      if nargin==1
         meas= double( vi(i).meas(idx) );
      else
         meas= double( vi(i).meas(idx) ) -  ...
               double( vh(1).meas(idx) );
      end
      vv(i) = eidors_obj('data',vi(i).name, ...
           'configuration', sprintf('Data from %dth stimulation',j), ...
           'meas', meas);
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
      meas= NaN*ones(ne*na,1); % Nan's will be replaced
      offset= (1:ne) - min(ne/2,i) + min(l_vi-i-ne/2, 0);
      for j= i+offset;
         f=rem(j-1,ne)+1; %extract jth frame
         idx= (f-1)*na+( 1:na );
         meas(idx) = vi(j).meas;
      end
      vv(i).meas= meas;
   end
