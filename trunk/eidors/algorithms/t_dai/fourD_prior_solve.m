function img= fourD_prior_solve( inv_model, data1, data2)
% fourD_prior_solve-- inverse solver to account for temporal
% and 3D spatial correlation
% img= fourD_prior_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2007, Tao Dai and Andy Adler. Licenced under the GPL Version 2
% $Id$

fwd_model= inv_model.fwd_model;
time_steps = inv_model.fourD_prior.time_steps;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'fourD_prior_solve');
if ~isempty(one_step_inv)
    eidors_msg('fourD_prior_solve: using cached value', 2);
else
    img_bkgnd = calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( fwd_model, img_bkgnd);

    %%%%%Calc temporal correlation%%%%%%%%%%%%%%%%%%%%%%%%%%
    [gama,Gama_x] = calc_Gama( inv_model, data1,data2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inv_model.fourD_prior.temporal_corr_mat = Gama_x;
    inv_model.fourD_prior.temporal_weight = gama;

    one_step_inv= data_form( inv_model, J );

    eidors_obj('set-cache', inv_model, 'fourD_prior_solve', one_step_inv);
    eidors_msg('fourD_prior_solve: setting cached value', 2);
end

if fwd_model.normalize_measurements
    dva= data2 ./ data1 - 1;
else
    dva= data2 - data1;
end
l_dva = size( dva, 2);

idx= [-time_steps:time_steps]'*ones(1,l_dva) + ...
    ones(time_steps*2+1,1)*(1:l_dva);
% replicate first and last measurements
idx(idx<1) = 1;
idx(idx>l_dva) = l_dva;

dvat= reshape(dva(:,idx),[],l_dva);
%%%%%% calc averaged measurement  %%%%%%%%%%%%%%%
% ts_vec = -time_steps:time_steps;
% weight = (gama.^abs(ts_vec));
% y_avg = dvat * weight(:) / sum(weight);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol = one_step_inv *dvat(:,1+time_steps); %y_avg;%

% create a data structure to return
img.name= 'solved by fourD_prior_solve';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% calculate the one_step_inverse using the data form
% CovX * J' * inv(J*CovX*J' + CovZ)
%   iRtR*Jt/(Ji*RtR*Jt +  hp^2*iW);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function one_step_inv= data_form( inv_model, J)
space_prior= inv_model.fourD_prior.space_prior;
time_weight= inv_model.fourD_prior.temporal_weight;
ts         = inv_model.fourD_prior.time_steps;

space_Reg= feval(space_prior, inv_model);

% %%%% calc space_Reg with inter-layer(z direction) coef considered

P = calc_P_3d(inv_model,inv(space_Reg));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iRtRJt_frac=  (P*J');
JiRtRJt_frac= J*iRtRJt_frac;

% JiRtRJt_mult accounts for different parts of the
% frame being taken at different times
delta_vec= calc_delta( inv_model, J);
delta_vec1= delta_vec*ones(1,length(delta_vec));
JiRtRJt_mult = time_weight.^abs(delta_vec1 - delta_vec1');

[x,y]= meshgrid(-ts:ts,  -ts:ts);
%     time_w_mat= time_weight.^abs(x-y)
time_w_mat= inv_model.fourD_prior.temporal_corr_mat;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate correlation matrix of images from measurements
function [gama,Gama] = calc_Gama( inv_model, data1,data2)
cov_n = inv_model.fourD_prior.cov_noise;%cov matrix of noise is
% calculated in advance in real meas,
% it can be obtained from baseline test.
ts = inv_model.fourD_prior.time_steps;

data = data2-data1;
Gama_y = corrcoef(data);%corr matrix of data
cov_y = cov(data');%cov matrix of data

JcovxJt = (cov_y-cov_n);%cov_y=JcovxJt+cov_n

% find out the best r to describe correlation between
% adjacent images. 
r=.01:.01:.99;
[k,l]= meshgrid(-ts:ts,-ts:ts);
I = eye(ts*2+1);
% cov_n_telda = kron(I,cov_n);
% cov_y_telda = kron(Gama_y,cov_y);%
MP = Gama_y*norm(cov_y,'fro') -I*norm( cov_n,'fro');%matrix precalculated
for i =1:length(r)
    G = r(i).^abs(k-l);
    e(i) = norm(MP - G*norm(JcovxJt,'fro'),'fro'); %this is equivalent to above. but save mem.
end
[e_min,ind] = min(e);
gama = r(ind)%correlation coefficient of adjacent img
Gama = gama.^abs(k-l);%correlation matrix of image set

function P_3d = calc_P_3d(inv_model,P_2d)

% P_C = exponential_covar_prior( inv_model );
P_C = calc_covar_prior( inv_model );%a simplified calculation of "exponential_covar_prior"

P_N = sqrt(P_2d);
P_3d = P_N*P_C*P_N;
