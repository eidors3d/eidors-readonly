function eidors_colourbar(max_scale,ref_lev, cb_shrink_move, greyscale)
% EIDORS_COLOURBAR - create an eidors colourbar with scaling to image
% usage: eidors_colourbar(max_scale,ref_lev)
%    ref_lev:   centre of the colour scale
%    max_scale: max difference from colour scale centre 
%
% Optional parameter:
%    cb_shrink_move(1) = horizontal shrink (relative)
%    cb_shrink_move(2) = vertial shrink (relative)
%    cb_shrink_move(3) = horizontal move (absolute screen units)
%
% The colorbars are removed with colorbar('delete')

% (C) 2005-2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

   hh= colorbar; 
   % make colourbar smaller and closer to axis
   if nargin >= 3
      posn= get(hh,'Position');
      cbsm = cb_shrink_move; 
      if ~all(cbsm == [1,1,0]); 
         posn = [posn(1) - cbsm(3), posn(2) + posn(4)*(1-cbsm(2))/2, ...
                 posn(3) * cbsm(1), posn(4) * cbsm(2)];
        
         set(hh,'Position', posn );

      end
   end

   %FIXME = AA+CG 30/1/12
   if nargin <4; 
       greyscale=[]; 
   else
       warning(['eidors_colourbar: greyscale is an experimental feature'...
           'and will be re-implemented']);
   end

   % Stop scale from being too small
   if max_scale<abs(ref_lev)
      if max_scale < 1e-10; max_scale = 1e-10; end
   else
      if max_scale/abs(ref_lev) < 1e-4; max_scale = ref_lev*1e-4; end 
   end

   % Get colormap limits  and move bottom so we don't see the background colour 
   ylim = get(hh,'Ylim');
   ylim(1)= ylim(1)+1;
   %FIXME = AA+CG 30/1/12
   if ~isempty(greyscale); ylim(1) = ylim(1) + 1 ; end
   set(hh,'Ylim',ylim);

   c_ctr = mean(ylim);
   c_max = ylim(2) - c_ctr;

%  in order to make the labels clean, we round to a near level
   OrdOfMag = 10^floor(log10(max_scale));
   scale_r  = OrdOfMag * floor( max_scale / OrdOfMag + 2*eps );
   ref_r = OrdOfMag * round( ref_lev / OrdOfMag );
   
   %FIXME = AA+CG 30/1/12
   
if isempty(greyscale)
    %   tick_vals = [-1:0.2:1]*max_scale + ref_r;
    tick_vals = [-1:0.5:1]*scale_r + ref_r;
else
%     tick_vals = [0:0.2:1]*max_scale;
     tick_vals = [0:0.2:1]*scale_r;
end

   % ref_lev goes to c_ctr. max_scale goes to c_max
%FIXME - need a switch to control use of max scale
   tick_locs = (tick_vals - ref_lev)/max_scale * c_max + c_ctr;
if isempty(greyscale) 
    tick_locs = (tick_vals - ref_lev)/max_scale * c_max + c_ctr;
else
    tick_locs = tick_vals*c_max*2/max_scale +2.5;
end

   set(hh,'YTick', tick_locs');
   set(hh,'YTickLabel', tick_vals');

