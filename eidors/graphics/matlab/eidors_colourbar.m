function hh= eidors_colourbar(max_scale,ref_lev, cb_shrink_move, greyscale)
% EIDORS_COLOURBAR - create an eidors colourbar with scaling to image
% usage: eidors_colourbar( img )
%    show a colourbar on the current axis representing img
%
% usage: hh= eidors_colourbar( img )
%    return a handle to colourbar created 
%
% usage: eidors_colourbar(max_scale,ref_lev)
%    ref_lev:   centre of the colour scale
%    max_scale: max difference from colour scale centre 
%
% Optional parameter:
%    cb_shrink_move(1) = horizontal shrink (relative)
%    cb_shrink_move(2) = vertial shrink (relative)
%    cb_shrink_move(3) = horizontal move (absolute screen units)
%
% Make sure you use 'axis tight' in order for cb_shrink_move to
%    give a correct looking image
%
% KNOWN ISSUE: if you use cb_shrink_move, then matlab will
%   forget the link between the figure and its colorbar. Future
%   plots in the same axis will continue to shrink. In general, the
%   axis will need to be cleared or reinitialized.
% EXAMPLE:
%   show_slices(img,2);
%    p = get(gca,'position') 
%   eidors_colourbar(img);
%    set(gca,'position',p); %%% Reset axes after colourbar and move
%
% The colorbars are removed with colorbar('delete')

% (C) 2005-2010 Andy Adler. License: GPL version 2 or version 3
% $Id$


if isstr(max_scale) && strcmp(max_scale,'UNIT_TEST'); do_unit_test; return; end

% if called as a simple eidors colourbar function
if isstruct(max_scale) && strcmp(max_scale.type,'image')
    calc_colours(max_scale,[],1);
    return
end

% Now deal with the other cases
   
   hh= colorbar; 
   % make colourbar smaller and closer to axis
   if nargin >= 3

      axpos = get(gca,'Position');
      posn= get(hh,'Position');
      cbsm = cb_shrink_move; 
      if ~all(cbsm == [1,1,0]); 
         posn = [posn(1) - cbsm(3), posn(2) + posn(4)*(1-cbsm(2))/2, ...
                 posn(3) * cbsm(1), posn(4) * cbsm(2)];
         set(hh,'Position', posn );
         set(gca,'Position',axpos);
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
   cmsz = size(colormap,1);
   cbsz = ylim(2) - ylim(1);
   unit = cbsz/cmsz;
   ylim(1)= ylim(1)+unit;
   %FIXME = AA+CG 30/1/12
   if ~isempty(greyscale); ylim(1) = ylim(1) + unit ; end
   set(hh,'Ylim',ylim);

   c_ctr = mean(ylim);
   c_max = ylim(2) - c_ctr;


% COMMENTS: AA - 4 apr 13
% If we want to ahve less aggressive rounding, we can do this
   F= 2;
   OrdOfMag = 10^floor(log10(max_scale*F))/F;

   scale1  = floor( max_scale / OrdOfMag + 2*eps );
% disp([scale1, OrdOfMag, max_scale]); => DEBUG
   if     (scale1/F >= 8);  fms = F*8;   ticks=[1]/2; 
   elseif (scale1/F >= 6);  fms = F*6;   ticks=[1]/2; 
   elseif (scale1/F >= 4);  fms = F*4;   ticks=[1]/2;
   elseif (scale1/F >= 3);  fms = F*3;   ticks=[1:2]/3;
   elseif (scale1/F >= 2);  fms = F*2;   ticks=[1]/2;
   elseif (scale1/F >= 1.5);fms = F*1.5; ticks=[1:2]/3;
   elseif (scale1/F >= 1);  fms = F*1;   ticks=[1]/2;
   else   (scale1/F >= 0.5);fms = F*0.5; ticks=[1]/2;
   end

   scale_r  = OrdOfMag * fms;

%  in order to make the labels clean, we round to a near level
   OrdOfMag = 10^floor(log10(max_scale));
   ref_r = OrdOfMag * round( ref_lev / OrdOfMag );
   
   %FIXME = AA+CG 30/1/12
if isempty(greyscale)
    %   tick_vals = [-1:0.2:1]*max_scale + ref_r;
    tick_vals = [-1,-fliplr(ticks),0,ticks,1]*scale_r + ref_r;
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

if fms==3; keyboard; end

   if nargin >= 3
% RESET OUR AXES
      if ~all(cbsm == [1,1,0]); 
         set(gca,'position',axpos);
      end
   end

   % Clear hh do it is not displayed, unless output requested
   if nargout ==0; clear hh; end

end

% DEBUG CODE ATTEMPTING TO FIX CB
function debug_code
         a = get(hh);
         set(hh,'Position', posn );
         a = rmfield(a,'CurrentPoint');
         a = rmfield(a,'TightInset');
         a = rmfield(a,'BeingDeleted');
         a = rmfield(a,'Type');
         a.Position = posn;
         set(hh,a);
         op = get(hh,'OuterPosition') 
         set(hh,'Position', posn );
         op1= get(hh,'OuterPosition') 
         set(hh,'OuterPosition',op);
         op2= get(hh,'OuterPosition') 
         set(gca,'Position',axpos);
end

function do_unit_test
   imgno = 1;
   imdl = mk_common_model('n3r2',[16,2]);
   img = mk_image(imdl);
   img=rmfield(img,'elem_data');
   img.node_data(1:252)= (0:251)/100 - 1.05;
   fprintf('CMAX = %4.2f CMIN = %4.2f\n', ...
       max(img.node_data), min(img.node_data) );


   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2); eidors_colourbar(img);

   % ref_level is not displayed, but will have its effect
   img.calc_colours.ref_level = -0.1234;
   img.calc_colours.ref_level =  0;
   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2); eidors_colourbar(img);

   img.calc_colours.cmax = 1;
   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2); eidors_colourbar(img);

   MV =-0.05;
   img.calc_colours.cb_shrink_move = [.5,.8,MV];
   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2);
   eidors_colourbar(img);

   subplot(3,3,imgno); imgno=imgno+1;
   show_fem(img,1);

   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2);
    p = get(gca,'position');
   eidors_colourbar(img);
    set(gca,'position',p);

   show_slices(img,2);
    eidors_colourbar(img);

   subplot(3,3,imgno); imgno=imgno+1;
   img.node_data = 26e1*abs(img.node_data);
   show_slices(img,2);
   eidors_colourbar(img);

   subplot(3,3,imgno); imgno=imgno+1;
   img.node_data(1:252)= abs( (0:251)/100 - 1.05 );
   show_slices(img,2);
   eidors_colourbar(img);


   subplot(3,3,imgno); imgno=imgno+1;
   img.calc_colours.ref_level = 0.5;
   img.calc_colours.clim      = 0.5;
   show_slices(img,2);
   eidors_colourbar(img);

end
