function Recon_GREIT_ndiff( savename );
% Reconstruct GREIT images using GREIT algorithm
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   [RM] = calc_RM( 5);
   normalize_flag = 1;
   save(savename, 'RM','normalize_flag');

function RM = calc_RM( noiselev)
   load sim_targets.mat

   D = desired_soln( xyzr_pt, 0.25 );
   Y = vi./(vh*ones(1,size(vi,2))) - 1; 
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
               (dx+x(:)+xyc(2,i)).^2 + (dy+y(:)+xyc(1,i)).^2 ...
                        < radius^2 );
         end
      end
%     PSF(:,i) = PSF(:,i)/sum(PSF(:,i));
   end
