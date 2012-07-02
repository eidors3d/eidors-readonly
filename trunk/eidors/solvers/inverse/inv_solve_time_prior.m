function img= inv_solve_time_prior( inv_model, data1, data2)
% INV_SOLVE_TIME_PRIOR inverse solver to account for time differences
% img= inv_solve_time_prior( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix
%
% Parameters
%   inv_model.inv_solve_time_prior.time_steps => time_steps

% TODO: This function really should be calling the proper
%   prior calculator functions, and not reimplementing
%   them internally

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

fwd_model= inv_model.fwd_model;
time_steps = inv_model.inv_solve_time_prior.time_steps;
l_ts  = time_steps*2 + 1;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'inv_solve_time_prior');
if ~isempty(one_step_inv)
    eidors_msg('inv_solve_time_prior: using cached value', 2);
else
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( fwd_model, img_bkgnd);

%   one_step_inv= standard_form( inv_model, J );
    one_step_inv= data_form( inv_model, J );

    eidors_obj('set-cache', inv_model, 'inv_solve_time_prior', one_step_inv);
    eidors_msg('inv_solve_time_prior: setting cached value', 2);
end

dva = calc_difference_data( data1, data2, inv_model.fwd_model);

l_dva = size( dva, 2);

idx= [-time_steps:time_steps]'*ones(1,l_dva) + ...
     ones(l_ts,1)*(1:l_dva);
% replicate first and last measurements
idx(idx<1) = 1;
idx(idx>l_dva) = l_dva;

dvat= reshape(dva(:,idx),[],l_dva);
 
sol = one_step_inv * dvat;

% create a data structure to return
img.name= 'solved by inv_solve_time_prior';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% calculate the one_step_inverse using the standard
% formulation (JtWJ + hp^2*RtR)\JtW
function one_step_inv= standard_form( inv_model, J )
    RtR = calc_RtR_prior( inv_model );
    W   = calc_meas_icov( inv_model );
    hp  = calc_hyperparameter( inv_model );

    time_steps = inv_model.inv_solve_time_prior.time_steps;
    l_ts  = time_steps*2 + 1;

    JtWJ = kron( speye(l_ts), J'*W*J);
    JtW  = kron( speye(l_ts), J'*W);
    one_step_inv= (JtWJ +  hp^2*RtR)\JtW;

    n_el= size(J,2);
    one_step_inv= one_step_inv(n_el*time_steps + (1:n_el),:);

% calculate the one_step_inverse using the data form
% CovX * J' * inv(J*CovX*J' + CovZ)
%   iRtR*Jt/(Ji*RtR*Jt +  hp^2*iW);
function one_step_inv= data_form( inv_model, J );
    space_prior= inv_model.prior_time_smooth.space_prior;
    time_weight= inv_model.prior_time_smooth.time_weight;
    ts         = inv_model.inv_solve_time_prior.time_steps;

    space_Reg= feval(space_prior, inv_model);

    iRtRJt_frac=  (space_Reg\J');
    JiRtRJt_frac= J*iRtRJt_frac;

    % JiRtRJt_mult accounts for different parts of the
    % frame being taken at different times
    delta_vec= calc_delta( inv_model, J);
    delta_vec1= delta_vec*ones(1,length(delta_vec));
    JiRtRJt_mult = time_weight.^abs(delta_vec1 - delta_vec1');

    [x,y]= meshgrid(-ts:ts,  -ts:ts);
    time_w_mat= time_weight.^abs(x-y);

    JiRtRJt= kron( time_w_mat, JiRtRJt_frac .* JiRtRJt_mult );

    %FIXME: do we multiply be JiRtRJt_mult here?
    iRtRJt=  kron( time_w_mat(ts+1,:), iRtRJt_frac );

    iW   = kron( speye(1+2*ts), inv( ...
                 calc_meas_icov( inv_model ) ));
    hp   = calc_hyperparameter( inv_model );

    one_step_inv= iRtRJt/(JiRtRJt +  hp^2*iW);

% if measurements are taken at different times,
% then calculate a delta of each wrt the centre time
function delta_vec= calc_delta( inv_model, J)
   stimulation= inv_model.fwd_model.stimulation;
   n_N= size(J,1);

   if isfield(stimulation(1),'delta_time')
      delta_time= [stimulation(:).delta_time];
      if diff(delta_time) ~= 0;
         error('All time steps must be same for kalman filter');
      end
   else
      delta_time=0;
   end

   % sequence is a vector location of each stimulation in the frame
   if delta_time == 0
      seq= size(J,1);
   else
      for i=1:length(stimulation)
         seq(i) = size(stimulation(i).meas_pattern,1);
      end
      seq= cumsum( seq );
   end

   delta_time= cumsum(delta_time);

   delta_vec= zeros(size(J,1),1);
   seq= [0;seq(:)];
   for i=1:length(seq)-1
      delta_vec( (seq(i)+1):seq(i+1) )= delta_time(i);
   end

   % normalize so middle time is centre, and max time is 1
   delta_vec= (delta_vec - mean(delta_vec)) / ...
              (sum(delta_time) + eps );

   
function do_unit_test
   time_steps=  3;
   time_weight= .8;

   [vh,vi,xyr_pt]=simulate_2d_movement(50);

   imdl_TS = mk_common_model( 'c2c2', 16 ); % 576 element
   imdl_TS.fwd_model.normalize_measurements= 0;
   imdl_TS.hyperparameter.value= 0.10;

   imdl_TS.RtR_prior= @prior_time_smooth;
   imdl_TS.prior_time_smooth.space_prior= @prior_noser;
   imdl_TS.prior_noser.exponent= .5;
   imdl_TS.prior_time_smooth.time_weight= time_weight;
   imdl_TS.prior_time_smooth.time_steps=  time_steps;
   imdl_TS.solve= @inv_solve_time_prior;
   imdl_TS.inv_solve_time_prior.time_steps=   time_steps;

image_select= length(xyr_pt)/2+1;; % this image is at 9 O'Clock
time_steps=  3; ts_expand= 5;
time_weight= .8;
ts_vec= -time_steps:time_steps;

    im_sel= image_select+ ts_vec*ts_expand;
    vi_sel= vi(:,im_sel);
    sel  = 1 + time_steps; % choose the middle

   img= inv_solve( imdl_TS, vh, vi_sel);
   img.elem_data= img.elem_data(:,sel);
   show_fem(img);




