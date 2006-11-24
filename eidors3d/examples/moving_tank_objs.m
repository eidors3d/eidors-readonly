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
% $Id: moving_tank_objs.m,v 1.18 2006-11-24 00:36:15 aadler Exp $

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

    case 22
        load iirc_data_2006
        vh= v_reference;
        vi= v_rotate(:,1:30);
        filename= 'iirc-around-clean';

    case 22.1
        load iirc_data_2006
        vh= v_reference;
        vi= v_rotate(:,1:10);
        filename= 'iirc-around-clean';

    case 23
        load iirc_data_2006
        vh= v_reference;
        vi= v_rotate(:,1:30);
        snr= norm(vi)/norm(vi - vh(:,ones(1,30)));
        rand('seed',50);
        vh= vh+snr*2*randn(size(vh));
        vi= vi+snr*2*randn(size(vi));
        filename= 'iirc-around-noise';

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
%          vi(i).meas = vi(i).meas + 25e-6*randn(size(vi(i).meas));
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

    case 2.1
        imdl= mk_common_model('c2c',16);
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
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';

        vi = extract_subframes(vi,4);
        vi = repack_subframes(vi,4);

    case 11.1
        imdl= mk_common_model('c2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= 'laplace_image_prior';
        imdl.solve= 'np_inv_solve';

        vi = extract_subframes(vi,4);
        vi = repack_subframes(vi,4);

    case 12
        imdl= mk_common_model('b2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= @laplace_image_prior;
        imdl.solve= @inv_kalman_diff;
        [imdl.fwd_model.stimulation(:).delta_time]=deal(1);

        vi = extract_subframes(vi,4);

    case 12.1
        imdl= mk_common_model('c2c',16);
        imdl.hyperparameter.value= 1e-2;
        imdl.RtR_prior= @laplace_image_prior;
        imdl.solve= @inv_kalman_diff;
        for i=1:16;
           imdl.inv_kalman_diff.sequence(i).meas_no= (i-1)*13+(1:13);
        end
        imdl.fwd_model = rmfield(imdl.fwd_model,'meas_select');
        imdl.fwd_model.normalize=0;

    case 14.1
        vi = sub_frame(vi,vh);
        vh = zeros(length(vi(1).meas),1);

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

        for i=1:16;
           imdl.time_prior_solve.sequence(i).meas_no= (i-1)*13+(1:13);
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

t=cputime;
imgs= inv_solve(imdl,vi,vh);
fprintf('solve time=%f\n',cputime-t);
animate_reconstructions(filename, imgs, clim);

% original = [ 1.1 2.1 3.1 4.1 5.1
%              1.2 2.2 3.2 4.2 5.2
%              1.3 2.3 3.3 4.3 5.3
%              1.4 2.4 3.4 4.4 5.4
%              1.5 2.5 3.5 4.5 5.5  ] 
% subseq = 2
% output   = [ 1.1 3.1 
%              1.2 4.2 
%              2.3 4.3 
%              2.4 5.4 
%              3.5 5.5 ] 
%
% Hardcoded for 16 electrode 2D adjacent stimulation
function ve= extract_subframes( vv, subseq)
   if isstruct(vv);
      vv= [vv(:).meas];
   end
   [st, els]= mk_stim_patterns(16, 1, '{ad}','{ad}', {}, 10);

   % thow away measurements at current elecs if needed
   if size(vv,1)==size(els,1)
      vv= vv(els,:);
   end

   ve= zeros( size(vv,1), floor(size(vv,2)*subseq/16) );
   dst=0; src=0;
   for k=0:16*size(ve,2)-1;
      pat = (1:13) + 13*rem(k,16);
      if 0==rem(k,16);     dst= dst+1; end
      if 0==rem(k,subseq); src= src+1; end
      ve(pat,dst) = vv(pat, src);
   end


% original = [ 1.1 2.1 3.1 4.1 5.1
%              1.2 2.2 3.2 4.2 5.2
%              1.3 2.3 3.3 4.3 5.3
%              1.4 2.4 3.4 4.4 5.4
%              1.5 2.5 3.5 4.5 5.5  ] 
% subseq = 2
%
% output   = [ 1.1 2.1 2.1 3.1 3.1 
%              1.2 2.2 2.2 2.2 3.2 
%              1.3 1.3 2.3 2.3 3.3 
%              1.4 1.4 2.4 2.4 2.4 
%              1.5 1.5 1.5 3.5 3.5 ] 
%
% Hardcoded for 16 electrode 2D adjacent stimulation

function ve = repack_subframes(vv,subseq);
   % duplicate first and last 
   ve= zeros( size(vv,1), floor(size(vv,2)*16/subseq) );
   vv= vv(:,[1:end,end]);
   dst=0; src=0;
   pat = (1:208)';
   for i=1:size(ve,2)
      for k=rem(subseq*i + (-subseq:-1),16);
         pat(13*k + (1:13)) = pat(13*k + (1:13)) + 208;
      end
      ve(:,i) = vv(pat);
   end

