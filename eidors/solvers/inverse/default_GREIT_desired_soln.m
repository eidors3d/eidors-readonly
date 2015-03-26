function PSF= default_GREIT_desired_soln(xyc, radius, opt)
% DEFAULT_GREI_DESIRED_SOLN The original desired solution for GREIT
% See also: CALC_GREIT_RM

% (C) 2009 Andy Adler. Licenced under GPL v2 or v3
% $Id$
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
   
end

