function eidors_colourbar(max_scale,ref_lev, cb_shrink_move)
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
   if nargin == 3
      posn= get(hh,'Position');
      cbsm = cb_shrink_move; 
      if ~all(cbsm == [1,1,0]); 
         posn = [posn(1) - cbsm(3), posn(2) + posn(4)*(1-cbsm(2))/2, ...
                 posn(3) * cbsm(1), posn(4) * cbsm(2)];
        
         set(hh,'Position', posn );

      end
   end

   % Get colormap limits  and move bottom so we don't see the background colour 
   ylim = get(hh,'Ylim');
   ylim(1)= ylim(1)+1;
   set(hh,'Ylim',ylim);

   c_ctr = mean(ylim);
   c_max = ylim(2) - c_ctr;

%  in order to make the labels clean, we round to a near level
   OrdOfMag = 10^floor(log10(max_scale));
   scale_r  = OrdOfMag * floor( max_scale / OrdOfMag );
   ref_r = OrdOfMag * round( ref_lev / OrdOfMag );

   tick_vals = [-1:0.5:1]*scale_r + ref_r;
   % ref_lev goes to c_ctr. max_scale goes to c_max
   tick_locs = (tick_vals - ref_lev)/max_scale * c_max + c_ctr;
   set(hh,'YTick', tick_locs');
   set(hh,'YTickLabel', tick_vals');

