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
   if ~isfield(opt, 'desired_solution_fn')
      D = default_GREIT_desired_soln( xyc, radius, opt);
   else
      D = feval(opt.desired_solution_fn, xyc, radius, opt);
   end

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
   RM = PJt/M;
   % This implements RM = D*Y'/(Y*Y');
   if 0
      Y = [Y, noiselev*eye(N_meas)];
      D = [D,          zeros(size(D,1),N_meas)];

      RMold = D*Y'/(Y*Y');
      if norm(RM-RMold,'fro')/norm(RM,'fro') > 1e-10; warning('not OK'); end
   end

function PSF= default_GREIT_desired_soln(xyc, radius, opt)
   c_obj = {xyc, radius, opt};
   PSF = eidors_obj('get-cache', c_obj, 'desired_solution');
   if ~isempty(PSF)
       return
   end
        
   xsz = opt.imgsz(1); ysz = opt.imgsz(2);
   sz= xsz * ysz;
   xmin = opt.meshsz(1); xmax = opt.meshsz(2);
   ymin = opt.meshsz(3); ymax = opt.meshsz(4);
   % scale radius to half the greater dimension
   radius = radius * 0.5 * max(xmax-xmin, ymax-ymin);
   [x,y]= ndgrid(linspace(xmin,xmax,xsz), linspace(ymin,ymax,ysz));
   x_spc = (xmax-xmin)/(xsz-1) * 0.5;
   y_spc = (ymax-ymin)/(ysz-1) * 0.5;
   PSF = zeros(sz,size(xyc,2));
   for i=1:size(xyc,2);
      for dx = linspace(-x_spc, x_spc, 5)
         for dy = linspace(-y_spc, y_spc, 5)
            PSF(:,i) = PSF(:,i) +  1/25*( ...
               (dx+x(:)-xyc(1,i)).^2 + (dy+y(:)-xyc(2,i)).^2 ...
                        < radius^2 );
         end
      end
%     PSF(:,i) = PSF(:,i)/sum(PSF(:,i));
   end
   eidors_obj('set-cache', c_obj, 'desired_solution', PSF);

   function opt = parse_options(opt)
       if ~isfield(opt, 'normalize'), opt.normalize = 1; end
       if ~isfield(opt, 'meshsz'),    opt.meshsz = [-1 1 -1 1]; end
       if ~isfield(opt, 'imgsz'),     opt.imgsz = [32 32]; end
       if ~isfield(opt, 'noise_covar'),
                     opt.noise_covar = 1;
       end
%   options.data_covariance [optional]
