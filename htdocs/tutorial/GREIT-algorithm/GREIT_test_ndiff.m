function [img,map]= GREIT_test_ndiff( ref_meas, reconst_meas )
% Reconstruct GREIT images using GREIT algorithm
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   [RM] = calc_RM(1, 12);

   % Expand ref_meas to the full size of reconst_meas
   num_meas = size(reconst_meas,2);
   ref_meas = ref_meas * ones(1,num_meas);
   dv = ( reconst_meas - ref_meas ) ./ ref_meas; % CHANGE IS HERE:

   % reconst image
   ds = RM*dv;

   img= reshape(ds, 32,32,num_meas);

function RM = calc_RM(data_file, noiselev)
   if data_file==1
      load sim_radmove_homog.mat
   else
      error('data_file not recognized');
   end

   D = desired_soln( xyzr_pt, 0.20 );
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
