function eidors_colourbar(max_scale,ref_lev)
% eidors_colourbar

% (C) 2005-2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

   hh= colorbar; delete(hh); hh=colorbar;
   % make colourbar smaller and closer to axis
   if sscanf(version,'%f',1) < 7.0
      p= get(hh,'Position');
      pm= p(2) + p(4)/2;
      set(hh,'Position', [p(1)+1.2*p(3), pm-p(4)*.6/2, p(3)*.6, p(4)*.6]);
   end

   % set scaling
%  lcm= size(colormap,1)/2+.5; - you would expect it to be this
   lcm= max(get(hh,'Ylim'))/2 + .5;
   OrdOfMag = 10^floor(log10(max_scale));
%  in order to make the labels clean, we round to a near level
   scale_r  = OrdOfMag * floor( max_scale / OrdOfMag );
%  ticks = lcm + (lcm-1)*[-1,0,+1]*scale_r/max_scale;
   ref_r = OrdOfMag * round( ref_lev / OrdOfMag );
   ofs      = [-1,0,+1]*scale_r + (ref_r-ref_lev);
   ticks = lcm + (lcm-1)*ofs/max_scale;
   set(hh,'YTick', ticks');
   set(hh,'YTickLabel', [-scale_r, 0, scale_r]'+ ref_r);

