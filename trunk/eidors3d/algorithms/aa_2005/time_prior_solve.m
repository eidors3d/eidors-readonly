function img= time_prior_solve( inv_model, data1, data2)
% TIME_PRIOR_SOLVE inverse solver to account for time differences
% img= aa_inv_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: time_prior_solve.m,v 1.4 2006-11-07 13:34:55 aadler Exp $

fwd_model= inv_model.fwd_model;
time_steps = inv_model.time_prior_solve.time_steps;
l_ts  = time_steps*2 + 1;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'time_prior_solve');
if ~isempty(one_step_inv)
    eidors_msg('time_prior_solve: using cached value', 2);
else
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( fwd_model, img_bkgnd);

%   one_step_inv= standard_form( inv_model, J );
    one_step_inv= data_form( inv_model, J );

    eidors_obj('set-cache', inv_model, 'time_prior_solve', one_step_inv);
    eidors_msg('time_prior_solve: setting cached value', 2);
end

% TODO: account for normalized data
if fwd_model.normalize_measurements
   dva= 1 - data2 ./ data1;
else   
   dva= data1 - data2;
end
l_dva = size( dva, 2);

idx= [-time_steps:time_steps]'*ones(1,l_dva) + ...
     ones(l_ts,1)*(1:l_dva);
% replicate first and last measurements
idx(idx<1) = 1;
idx(idx>l_dva) = l_dva;

dvat= reshape(dva(:,idx),[],l_dva);
 
sol = one_step_inv * dvat;

% create a data structure to return
img.name= 'solved by time_prior_solve';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% calculate the one_step_inverse using the standard
% formulation (JtWJ + hp^2*RtR)\JtW
function one_step_inv= standard_form( inv_model, J )
    RtR = calc_RtR_prior( inv_model );
    W   = calc_meas_icov( inv_model );
    hp  = calc_hyperparameter( inv_model );

    time_steps = inv_model.time_prior_solve.time_steps;
    l_ts  = time_steps*2 + 1;

    JtWJ = kron( speye(l_ts), J'*W*J);
    JtW  = kron( speye(l_ts), J'*W);
    one_step_inv= (JtWJ +  hp^2*RtR)\JtW;

    n_el= size(J,2);
    one_step_inv= one_step_inv(n_el*time_steps + (1:n_el),:);

% calculate the one_step_inverse using the standard
% formulation (JtWJ + hp^2*RtR)\JtW
function one_step_inv= data_form( inv_model, J );
    space_prior= inv_model.time_smooth_prior.space_prior;
    time_weight= inv_model.time_smooth_prior.time_weight;
    ts         = inv_model.time_prior_solve.time_steps;

    if isfield(inv_model.time_prior_solve,'inv_kalman_diff')
       sequence= inv_model.time_prior_solve.sequence;
    else
       sequence.meas_no= 1:size(J,1); % all data all the time
    end

    space_Reg= feval(space_prior, inv_model);

    iRtRJt_frac=  (space_Reg\J');
    JiRtRJt_frac= J*iRtRJt_frac;

    [x,y]= meshgrid(-ts:ts,  -ts:ts);
    time_w_mat= time_weight.^abs(x-y);

    JiRtRJt= kron( time_w_mat, JiRtRJt_frac );
    iRtRJt=  kron( time_w_mat(ts+1,:), iRtRJt_frac );

    iW   = kron( speye(1+2*ts), inv( ...
                 calc_meas_icov( inv_model ) ));
    hp   = calc_hyperparameter( inv_model );

    one_step_inv= iRtRJt/(JiRtRJt +  hp^2*iW);
