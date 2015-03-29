function [RM, PJt, M] = calc_GREIT_RM(vh,vi, xyc, radius, weight, options)
% CALCULATE GREIT reconstruction matrix
%   RM= calc_GREIT_RM(vh,vi, xyc, radius, weight, normalize)
% 
% Input:
%   vh     = homogeneous (reference) training measurements 
%   vi     = inhomogeneous training measurements 
%   xyc    = x,y position of targets 2xN
%   radius = requested resolution matrix
%      if scalar - apply resolution to all images:  recommend 0.25 for 16 elecs
%   weight = weighting matrix
%      if scalar   = weighting of noise vs signal
%      if 32^2 x N = weighting of each image output
%   options.normalize = 
%            0 -> regular difference EIT
%            1 -> normalized difference EIT
%   options.meshsz = [xmin xmax ymin ymax] size of mesh
%      defaults to [-1 1 -1 1]
%   options.imgsz  = [xsz ysz] size of the reconstructed image in pixels
%   options.noise_covar [optional]
%      covariance matrix of data noise
%   options.desired_solution_fn
%      specify a function to calculate the desired image. 
%      It must have the signature:
%      D = my_function( xyc, radius, options);
%      uses eidors_defualt('get','GREIT_desired_img') if not specified
% 
%
% (C) 2009 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   if ~isstruct(options)
       options.normalize = options;
   end
   opt = parse_options(options);

   if opt.normalize
      Y = vi./(vh*ones(1,size(vi,2))) - 1; 
   else
      Y = vi - (vh*ones(1,size(vi,2)));
   end
   try 
       f = opt.desired_solution_fn;
   catch
       f = eidors_default('get','GREIT_desired_img');
   end

   D = feval(f, xyc, radius, opt);
   
   if size(weight)==[1,1] % Can't use isscalar for compatibility with M6.5
       [RM, PJt, M] = calc_RM(Y,D,weight, opt);
   else
       error('not coded yet');
   end

function [RM, PJt, M] = calc_RM(Y, D, noiselev, opt)

   noiselev = noiselev * mean(abs(Y(:)));
   % Desired soln for noise is 0
   N_meas = size(Y,1);

   % This implements RM = D*Y'/(J*Sx*J + Sn);
   Sn = speye(N_meas) .* opt.noise_covar; % Noise covariance
   PJt= D*Y';
   M  = (Y*Y' + noiselev^2*Sn);
   RM =  left_divide(M',PJt')';    %PJt/M;
   % This implements RM = D*Y'/(Y*Y');
   if 0
      Y = [Y, noiselev*eye(N_meas)];
      D = [D,          zeros(size(D,1),N_meas)];

      RMold = D*Y'/(Y*Y');
      if norm(RM-RMold,'fro')/norm(RM,'fro') > 1e-10; warning('not OK'); end
   end



   function opt = parse_options(opt)
       if ~isfield(opt, 'normalize'), opt.normalize = 1; end
       if ~isfield(opt, 'meshsz'),    opt.meshsz = [-1 1 -1 1]; end
       if ~isfield(opt, 'imgsz'),     opt.imgsz = [32 32]; end
       if ~isfield(opt, 'noise_covar'),
                     opt.noise_covar = 1;
       end
%   options.data_covariance [optional]
