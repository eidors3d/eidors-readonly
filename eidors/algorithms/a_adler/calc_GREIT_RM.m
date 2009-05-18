function RM= calc_GREIT_RM(vh,vi, xyc, radius, weight, normalize)
% CALCULATE GREIT reconstruction matrix
%   RM= calc_GREIT_RM(vh,vi, xyc, radius weight, normalize)
% 
% Input:
%   vh     = homogeneous (reference) training measurements 
%   vi     = inhomogeneous training measurements 
%   xyc    = x,y position of targets 2xN
%   radius = requested weighting matrix  (recommend 0.25 for 16 electrodes)
%   weight = weighting matrix
%      if scalar   = weighting of noise vs signal
%      if 32^2 x N = weighting of each image output
%   normalize = 0 -> regular difference EIT
%                 -> normalized difference EIT
% 
%
% (C) 2009 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   if normalize
      Y = vi./(vh*ones(1,size(vi,2))) - 1; 
   else
      Y = vi - (vh*ones(1,size(vi,2)));
   end

   D = desired_soln( xyc, radius );

   if size(weight)==[1,1] % Can't use isscalar for compatibility with M6.5
       [RM] = calc_RM(Y,D,weight);
   else
       error('not coded yet');
   end

function RM = calc_RM(Y, D, noiselev)

   noiselev = noiselev * mean(abs(Y(:)));
   % Desired soln for noise is 0
   Y = [Y, noiselev*eye(208)];
   D = [D,          zeros(32^2,208)];

   RM = D*Y'/(Y*Y');

function PSF= desired_soln(xyc, radius)
   xsz= 32; ysz= 32; sz= xsz * ysz;
   lim= 1.00;
   [x,y]= ndgrid(linspace(-lim,lim,xsz), linspace(-lim,lim,ysz));
   spc = 2*lim/(xsz-1) * 0.5;
   PSF = zeros(sz,size(xyc,2));
   for i=1:size(xyc,2);
      for dx = linspace(-spc, spc, 5)
         for dy = linspace(-spc, spc, 5)
            PSF(:,i) = PSF(:,i) +  1/25*( ...
               (dx+x(:)-xyc(1,i)).^2 + (dy+y(:)-xyc(2,i)).^2 ...
                        < radius^2 );
         end
      end
%     PSF(:,i) = PSF(:,i)/sum(PSF(:,i));
   end

