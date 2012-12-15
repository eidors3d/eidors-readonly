function [NF,SE] = calc_noise_figure( inv_model, hp, iterations)
% CALC_NOISE_FIGURE: calculate the noise amplification (NF) of an algorithm
% [NF,SE] = calc_noise_figure( inv_model, hp, iterations)
%    inv_model  => inverse model object
%    hp         => value of hyperparameter to use (if not spec
%         then use the value of inv_model.hyperparameter.value)
%    iterations => number of iterations (default 10)
%       (for calculation of noise figure using random noise)
%    if hp is [] => use the hp on the inv_model
%  NF = calculated NF. SE = standard error on NF
%
% hp is specified, it will be used for the hyperparameter.
%    Otherwise the inv_model.hyperparameter will be used.
%
% Noise Figure must be defined for a specific measurement
% In order to specify data, use
%     inv_model.hyperparameter.tgt_data.meas_t1
%     inv_model.hyperparameter.tgt_data.meas_t2
%   to use a temporal solver (or the Kalman filter), the
%   measurement to perform the NF calc must also be specified,
%   using:
%     inv_model.hyperparameter.tgt_data.meas_select
%   otherwise, the middle measurement will be used
%
% In order to automatically simulate data, specify tgt_elems,
%   containing a vector of elements to use
% 
%     inv_model.hyperparameter.tgt_elems
%
% [NF,SE] = calc_noise_figure( inv_model, vh, vi)
%    interface provided for compatibility with the deprecated
%    calc_noise_params

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% A normal definition of noise power is based on power:
%      NF = SNR_z / SNR_x
%    SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
%    SNR_x = sumsq(x0) / var(x) = sum(x0.^2) / trace(ABRnB'A')
% but his doesn't work, because the signal spreads like
% the amplitude, not like the power, thus
%      NF = SNR_z / SNR_x
%    SNR_z = sum(|z0|/len_z) / std(z/len_z)
%    SNR_x = sum(|x0|/len_x) / std(x/len_x)

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'), do_unit_test, return, end

if nargin>=2 && numel(hp) == 1
   inv_model.hyperparameter.value= hp;
% Remove function parameter because it will recurse
   try; inv_model.hyperparameter = rmfield(inv_model.hyperparameter,'func'); end
end
if nargin == 3 && numel(hp) > 1
    h_data = hp;
    c_data = iterations;
else
    [inv_model, h_data, c_data] = process_parameters( inv_model );
end

%NF= nf_calc_use_matrix( inv_model, h_data, c_data);
%NF= nf_calc_iterate( inv_model, h_data, c_data); 
if nargin<3; iterations= 10; end
solver = inv_model.solve;
if isstr(solver) && strcmp(solver, 'eidors_default')
    solver = eidors_default('get','inv_solve');
end
if isa(solver,'function_handle')
    solver = func2str(solver);
end
switch solver
    case {'inv_solve_backproj'
            'inv_solve_conj_grad'
            'inv_solve_diff_GN_one_step'
            'inv_solve_trunc_iterative'
            'inv_solve_TSVD'
            'solve_use_matrix'}
        NF = nf_calc_linear( inv_model, h_data, c_data);
        SE = 0;
    otherwise
        [NF,SE]= nf_calc_random( inv_model, h_data, c_data, iterations);
end
eidors_msg('calculating NF=%f', NF, 2);

function [inv_model, h_data, c_data] = process_parameters( inv_model );

   if     isfield(inv_model.hyperparameter,'tgt_elems')
      [h_data, c_data]= simulate_targets( inv_model.fwd_model, ...
           inv_model.hyperparameter.tgt_elems);
   elseif isfield(inv_model.hyperparameter,'tgt_data')
      tgt_data= inv_model.hyperparameter.tgt_data;
      h_data= tgt_data.meas_t1;
      c_data= tgt_data.meas_t2;
   else
      error('unsure how to get data to measure signal');
   end

   % if hp is specified, then use that value

   if isfield(inv_model.hyperparameter,'func')
      funcname= inv_model.hyperparameter.func;
      if strcmp( class(funcname), 'function_handle')
         funcname= func2str(funcname);
      end

      if strcmp(funcname, 'choose_noise_figure')
         error('specifying inv_model.hp.func = choose_noise_figure will recurse');
      end
   end

function params = nf_calc_linear(imdl, vh, vi )
% params = GREIT_noise_params(imdl, homg_voltage, sig_voltage)
%  params(1,:) = Noise Figure = SNR(image) / SNR(data)
%
%  see also: eval_GREIT_fig_merit or using test_performance

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% NOTE THAT WE ASSUME A LINEAR ALGORITHM FOR THIS MEASURE

if 0 % old code with random noise
     % Keep this to validate while we test it
   Nnoise = 1000;
   noise = 0.01*std(vh)*randn(size(vh,1),Nnoise);
   vhn= vh*ones(1,Nnoise) + noise;
else % use independent noise model on each channel
   noise = 0.01*std(vh)*speye(size(vh,1));
   vhn= vh*ones(1,size(vh,1)) + noise;
end

signal_y = calc_difference_data( vh, vi,  imdl.fwd_model);
noise_y  = calc_difference_data( vh, vhn, imdl.fwd_model);

signal_x = inv_solve(imdl, vh, vi);  signal_x = signal_x.elem_data;
noise_x  = inv_solve(imdl, vh, vhn); noise_x  = noise_x.elem_data;

use_rec = 1;
try 
   use_rec = ~imdl.prior_use_fwd_not_rec;
end
if use_rec
   try
      VOL = get_elem_volume(imdl.rec_model);
   catch
      VOL = get_elem_volume(imdl.fwd_model);
   end
else
   VOL = get_elem_volume(imdl.fwd_model);
end
VOL = spdiags(VOL,0, length(VOL), length(VOL));

signal_x = VOL*signal_x;
noise_x = VOL*noise_x;

signal_x = mean(abs(signal_x),1);
noise_x  = mean(std(noise_x));
snr_x = signal_x / noise_x;

signal_y = mean(abs(signal_y),1);
noise_y  = mean(std(noise_y)); 
snr_y = signal_y / noise_y;

params= [snr_y(:)./snr_x(:)]';

eidors_msg('NF= %f', params, 1);

% NOTES on the calculations: AA - Feb 20, 2012
% SNR = mean(abs(x)); VAR = 
% ym= E[y]                
% Sy= E[(y-ym)*(y-ym)'] = E[y*y'] - ym*ym'
% ny = sqrt(trace(Sy))
% xm= E[x]  = E[R*y] = R*E[y] = R*ym
% Sx= E[(x-xm)*(x-xm)'] = E[x*x'] - xm*xm'
%   = E[R*ym*ym'*R'] = R*E[ym*ym']*R' = R*Sy*R'
% nx = sqrt(trace(Sx))
% 
% signal = mean(abs(x));
% 
% In this case, these are exactly the same:
%    noise_x  = mean(std(noise_x));
%    noise_x  = sqrt(mean(noise_x.^2,2));
   
   
function NF= nf_calc_use_matrix( inv_model, h_data, c_data)
% To model std(z) we use z=z0+n
% so that std(z) = sqrt(var(z)) = sqrt(1/L * E[n'*n])
% we know a priori that the mean noise is zero, thus
% we do not need to divide by L-1 in the variance
% E[n'*n] = E[ trace(n*n') ] = trace( cov_N )
% Thus, std(z) = sqrt( trace( cov_N )/L ) 
%              = sqrt( mean( diag( cov_N )))
% And,  std(x) = sqrt( mean( diag( cov_X )))
%
% To model cov_N, we consider independant noise
%  on each channel, cov_N= N*N', N=diag( sigma_i )
% And,              cov_X= X*X', X=reconst(N)
% where X= reconst(mdl, z0,z0+N)
%
% To run efficiently mean(diag(cov_N))=mean(sum(N.^2,2))
% The sum over X is actually weighted by the volume of
%  each element, so VOL.^2*(sum(X.^2,2)
   try 
       VOL = get_elem_volume(inv_model.rec_model)';
   catch
       VOL = get_elem_volume(inv_model.fwd_model)';
   end

   % calculate signal
   d_len   = size(h_data,1);
   delta   = 1e-2* mean(h_data);
   c_noise = c_data*ones(1,d_len) + eye(d_len);
   h_full  = h_data*ones(1,d_len);

   sig_data = mean(abs( ...
         calc_difference_data( h_data, c_data , inv_model.fwd_model ) ...
                       ));
   var_data = mean(sum( ...
         calc_difference_data( h_full, c_noise, inv_model.fwd_model ) ...
                       .^2, 2)); 
                      

   % calculate image 
   % Note, this won't work if the algorithm output is not zero biased
   [img0, img0n] = get_images( inv_model, h_data, c_data, ...
                               h_full, c_noise);

   i_len = length(img0);
   sig_img= VOL*abs(img0) / i_len;;
   var_img= VOL.^2*sum(img0n.^2 ,2) / i_len;
   
   NF = ( sig_data/ sqrt(var_data) ) / ( sig_img / sqrt(var_img)  );

   % For the record, the expression for var_img is derived as:
   % Equiv expresssions for var_img % given: A= diag(pp.VOLUME);
   % var_img= trace(A*RM*diag(Rn.^2)*RM'*A');
   % vv=A*RM*diag(Rn);var_img=trace(vv*vv'); var_img= sum(sum(vv.^2));
   % var_img= VOL2* (RM*diag(Rn)).^2
   % var_img= VOL2* RM.^2 * Rn.^2

   
% simulate homg data and a small target in centre
function [h_data, c_data]= simulate_targets( fwd_model, ctr_elems)

   homg= 1; % homogeneous conductivity level is 1

   %Step 1: homogeneous image
   sigma= homg*ones( size(fwd_model.elems,1) ,1);

   img= eidors_obj('image', 'homogeneous image', ...
                   'elem_data', sigma, ...
                   'fwd_model', fwd_model );
   h_data=fwd_solve( img );
   h_data= h_data.meas;

   %Step 1: inhomogeneous image with contrast in centre
   delta = 1e-2;
   sigma(ctr_elems) = homg*(1 + delta);
   img.elem_data = sigma;
   c_data=fwd_solve( img );
   c_data= c_data.meas;

function [img0, img0n] = get_images( inv_model, h_data, c_data, ...
                               h_full, c_noise);
   if isa(inv_model.solve,'function_handle')
      solve= func2str(inv_model.solve);
   else
      solve= inv_model.solve;
   end

% Test for special functions and solve them specially
   switch solve
   case 'ab_tv_diff_solve'
      error('Dont know how to calculate TV noise figure')

   case 'inv_solve_diff_kalman'
      inv_model.inv_solve_diff_kalman.keep_K_k1= 1;
      stablize = 6;
      img0 = inv_solve( inv_model, h_data, ...
                                   c_data*ones(1,stablize) );
      K= img0.inv_solve_diff_kalman.K_k1;
      img0.elem_data = K*calc_difference_data( h_data , c_data , inv_model.fwd_model);
      img0n.elem_data= K*calc_difference_data( h_full , c_noise, inv_model.fwd_model);

   otherwise
      img0 = inv_solve( inv_model, h_data, c_data);
      if nargin>4
      img0n= inv_solve( inv_model, h_full, c_noise);
      end
   end

   % Need elem or nodal data
   if isfield(img0,'node_data');
      img0 = img0.node_data;
   else
      img0 = img0.elem_data;
   end

   if isfield(img0n,'node_data');
      img0n = img0n.node_data;
   else
      img0n = img0n.elem_data;
   end

% OLD CODE - iterate
function NF= nf_calc_iterate( inv_model, h_data, c_data);
   try 
       VOL = get_elem_volume(inv_model.rec_model)';
   catch
       VOL = get_elem_volume(inv_model.fwd_model)';
   end
   % calculate signal
   d_len   = size(h_data,1);
   delta   = 1e-2* mean(h_data);
   sig_data = mean(abs( ...
         calc_difference_data( h_data, c_data , inv_model.fwd_model ) ...
                       ));
   % calculate image 
   % Note, this won't work if the algorithm output is not zero biased

   [img0] = get_images( inv_model, h_data, c_data);
%  sig_img= mean(VOL'.*abs(img0.elem_data));
   sig_img= VOL*abs(img0) / length(img0);

   % Now do noise
   var_data= 0;
   var_img = 0;
   for i=1:d_len
      this_noise = -ones(d_len, size(c_data,2))/(d_len-1);
      this_noise(i,:) = 1;
      c_noise = c_data + this_noise;
      [imgn] = get_images( inv_model, h_data, c_noise);
      if 1
         var_data = var_data + mean(sum( ...
            calc_difference_data( h_data, c_noise, inv_model.fwd_model ) ...
                          .^2, 2)); 
%        var_img= var_img +  mean( (VOL'.*imgn.elem_data).^2 ); 
         var_img= var_img +  (VOL.^2)*sum(imgn.elem_data.^2,2 ) / length(imgn.elem_data); 
      else
         % OLD APPROACH BASED ON variance, rather than matrix calcs
         var_data = var_data + var( ...
            calc_difference_data( h_data, c_noise, inv_model.fwd_model ) ...
                                 ); 
         var_img= var_img + var( VOL'.*imgn.elem_data ); 
      end
   end
   var_data = var_data / d_len;
   var_img  = var_img  / d_len;
   NF = ( sig_data/ sqrt(var_data) ) / ( sig_img / sqrt(var_img)  );

function [NF,SE]= nf_calc_random( rec, vh, vi, N_RUNS);
   eidors_cache('boost_priority',-2); % low priority values

   imgr= inv_solve(rec, vh, vi);

   if isfield(imgr,'node_data');
      img0 = imgr.node_data;
      try
          VOL = get_elem_volume(rec.rec_model, 1);
      catch
          VOL = get_elem_volume(rec.fwd_model, 1);
      end
   else
      n_els = length(rec.fwd_model);
      img0 = imgr.elem_data(1:n_els); % elem_data can also contain movement
      try
          VOL = get_elem_volume(rec.rec_model, 0);
      catch
          VOL = get_elem_volume(rec.fwd_model, 0);
      end
   end

   sig_ampl = mean( abs( VOL .* img0 )) / ...
              mean( abs( calc_difference_data( vh, vi, rec.fwd_model )));

% Estimate Signal Amplitude
   for i=1:N_RUNS
      vn= addnoise(vh, vi, 1.0);

      imgr= inv_solve(rec, vh, vn);

      if isfield(imgr,'node_data'); img0 = imgr.node_data;
      else;                         img0 = imgr.elem_data(1:n_els);
      end

      noi_imag(i) = std( VOL .* img0 );
      noi_sgnl(i) = std( calc_difference_data( vh, vn, rec.fwd_model ));
   end

   noi_ampl = noi_imag./noi_sgnl;
   NF =  mean(noi_ampl/sig_ampl);
   SE =  std(noi_ampl/sig_ampl)/sqrt(N_RUNS);
   eidors_msg('NF= %f+/-%f', NF, SE, 1);

   eidors_cache('boost_priority',2);

function noise= addnoise( vh, vi, SNR);
      if isstruct(vh); vh= vh.meas; end
      if isstruct(vi); vi= vi.meas; end
      noise = randn(size(vh));
      noise = noise*std(vh-vi)/std(noise);
      noise = vh + SNR*noise;

function do_unit_test
    ll = eidors_msg('log_level',1);
    test1; % can we deal with c2f ?
    imdl = mk_common_model('a2t2',16); test2(imdl);
    imdl = mk_common_model('d2t2',16); test2(imdl);
    imdl = mk_common_model('a2c2',16); test2(imdl);
    imdl = mk_common_model('d2c2',16); test2(imdl);
    ll = eidors_msg('log_level',ll);

function test1
    % big model with c2f and supposedly an NF of 0.5
    fmdl = mk_library_model('pig_23kg_16el');
    [fmdl.stimulation fmdl.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}');
    fmdl = mdl_normalize(fmdl, 1);  % Use normalized difference imaging
    opt.noise_figure = 0.5; opt.imgsz = [64 64];
    imdl = mk_GREIT_model(fmdl, 0.25, [], opt);
    % homogeneous measurement
    img = mk_image(fmdl,1);
    vh = fwd_solve(img);
    % inhomogeneous measurement
    select_fcn = inline('(x-0).^2+(y-0).^2+(z-0.5).^2<0.1^2','x','y','z');
    mfrac = elem_select(fmdl, select_fcn);
    img.elem_data = img.elem_data + mfrac*0.1;
    vi = fwd_solve(img);
    
    nf1 = calc_noise_params(imdl, vh.meas, vi.meas);
    
    % We use different measurements so naturally we don't get 0.5 here
    imdl.hyperparameter.tgt_data.meas_t1 = vh.meas;
    imdl.hyperparameter.tgt_data.meas_t2 = vi.meas;
    try
        % calc_noise_figure doens't support dual models
        nf2 = calc_noise_figure(imdl);
    catch
        nf2 = 0;
    end
    unit_test_cmp('Noise fig implementations',nf1, nf2, 1e-2);

function test2(imdl)
    fmdl = imdl.fwd_model;
    % homogeneous measurement
    img = mk_image(fmdl,1);
    vh = fwd_solve(img);
    % inhomogeneous measurement
    select_fcn = inline('(x-0).^2+(y-0).^2.^2<15^2','x','y','z');
    mfrac = elem_select(fmdl, select_fcn);
    img.elem_data = img.elem_data + mfrac*0.1;
    vi = fwd_solve(img);

    nf1 = calc_noise_params(imdl, vh.meas, vi.meas);
    eidors_msg(nf1,0);

imdl.hyperparameter.tgt_data.meas_t1 = vh.meas;
imdl.hyperparameter.tgt_data.meas_t2 = vi.meas;
% calc_noise_figure doens't support dual models
nf2 = calc_noise_figure(imdl,[],1000);
unit_test_cmp('Noise fig implementations',nf1, nf2, 1e-2);

