function colours= calc_colours(img, scale, do_colourbar)
% colours= calc_colours(img, scale, do_colourbar)
% Calculate a colour for each image element 
% Conductive (positive) areas are shown in red
% Non-Conductive (negative) areas are shown in blue
%
% Usage:
%   c_img = calc_colours( r_img);
%   image( c_img );
% 
% img - an eidors image object, OR
%     - a 2D image matrix
%
% scale - colour value corresponding to maximum
%       - if not specified or scale==[] => autoscale
%
%    When autoscale is set, an appropriate background
%    reference conductivity is selected, if possible
%
% do_colourbar ==1 => show a Matlab colorbar with
%        appropriate scaling
%
% Colour maps are controlled by the global variable
%   eidors_colours. The following settings are defaults
%
%   eidors_colours.greylev = .2;
%      greylev is the colour of the ref_level (the non-changing regions
%      for difference imaging). Negative values indicate black (inversed
%      colour). For almost white, greylev=.01; Black=> greylev=-.01
%   eidors_colours.sat_adj = .9;
%       max G,B when R=1
%   eidors_colours.backgnd= [.5,.5,.15]; 
%      colour for non image regions, ie. the border around the image
%   eidors_colours.ref_level = 0  
%      conductivity of this value is centre of colour mapping. Normally,
%      this would be set to the conductivity of the background
%   eidors_colours.mapped_colour= 0; % use colormap function
%      if mapped_colour is non-zero, it indicates the colourmap
%      size; otherwise, RGB values are used.
%               Total colourmap is 2*mapped_colour
%      using mapped_colour allows matlab to print vector graphics to eps
%
% These global values may also be set via calc_colours
%    eg. calc_colours('ref_level',1)
%
%    eg. calc_colours('mapped_colour',256)
%            Use this to allow printing vector eps files to get around
%            a matlab bug with warning 'RGB color data not ...'
%
% Values may be obtained by calling calc_colours with 1 parameter
%   r_lev= calc_colours('ref_level');
%

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_colours.m,v 1.25 2006-07-29 15:50:26 aadler Exp $  

% Process input args
if nargin==0; return; end

% Now process scaling 
autoscale=1;
if nargin >= 2; if ~isempty(scale)
    autoscale=0;
end; end

% Do we want a colourbar
if nargin < 3;
   do_colourbar = 0;
end

pp=get_colours;
 
if isstr(img)
    global eidors_colours;
    if nargin==1
       colours = getfield(eidors_colours, img);
    else
       eidors_colours = setfield(eidors_colours, img, scale);
    end
    return;
elseif isfield(img,'type')
   if strcmp( img.type, 'image' )
      elem_data= img.elem_data(:); %col vector
   else
      error('calc_colours: input is not eidors image object');
   end
else
   elem_data= img(:);
end

if isempty(elem_data)
    colours = 'k'; %black
    return;
end

% remove background
e= length(elem_data);
if autoscale
    % we attempt to extimate the mode.
    % find the mean of the most common 50% of
    % elements
    s_ed= sort(elem_data);
    elem_data = elem_data - mean(s_ed( ceil(.25*e):floor(.75*e) ));
keyboard
else
    elem_data = elem_data - pp.ref_level;
end

backgnd= isnan(elem_data);
elem_data(backgnd)= mean( elem_data(~backgnd));

if autoscale
   scale =  max(abs(elem_data)) + eps;
end

if ~pp.mapped_colour
   [red,grn,blu] = blu_red_axis( pp, elem_data / scale, backgnd );
   colours= shiftdim( [red,grn,blu], -1);
else
   colours=set_mapped_colour(pp, scale, backgnd, elem_data);
end

% print colorbar if do_colourbar is specified
if do_colourbar
   if ~pp.mapped_colour
       warning('Colorbar not available without mapped_colour option');
   else

   hh= colorbar;
   p= get(hh,'Position');
   pm= p(2) + p(4)/2;
   set(hh,'Position', [p(1)+p(3), pm-p(4)*.6/2, p(3)*.6, p(4)*.6]);

   % set scaling
   lcm= size(colormap,1)/2+.5;
   OrdOfMag = 10^floor(log10(scale));
   scale_r  = OrdOfMag * floor( scale / OrdOfMag );
   ticks = lcm + (lcm-1)*[-1,0,+1]*scale_r/scale;
   set(hh,'YTick', ticks');
   set(hh,'YTickLabel', [-scale_r, 0, scale_r]'+ pp.ref_level);

   end
end


%scaled data must go from -1 to 1
function [red,grn,blu] = blu_red_axis( pp, scale_data, backgnd )
   D= sign(pp.greylev+eps); %force 0 to 1
   glev= abs(pp.greylev);
   F= 3*pp.sat_adj;

   red= D*F*abs(scale_data+D/F) - D + (D==-1);
   red= red.*(red>0).*(red<1) + (red>=1);
   red= red*(1-glev) + glev;

   grn= D*F*abs(scale_data    ) - D + (D==-1);
   grn= grn.*(grn>0).*(grn<1) + (grn>=1);
   grn= grn*(1-glev) + glev;

   blu= D*F*abs(scale_data-D/F) - D + (D==-1);
   blu= blu.*(blu>0).*(blu<1) + (blu>=1);
   blu= blu*(1-glev) + glev;

   red(backgnd) = pp.backgnd(1);
   grn(backgnd) = pp.backgnd(2);
   blu(backgnd) = pp.backgnd(3);


function pp=get_colours;
   global eidors_colours;

   if ~isfield( eidors_colours, 'greylev' );
      eidors_colours.greylev = -.01;
   end
   if ~isfield( eidors_colours, 'sat_adj' );
      eidors_colours.sat_adj = .9;
   end
   if ~isfield( eidors_colours, 'backgnd' );
      eidors_colours.backgnd= [.5,.5,.15];
   end
   % better to set the default to mapped_colour. Matlab
   %  seems to like this better anyway (ie less bugs)
   if ~isfield( eidors_colours, 'mapped_colour' );
      eidors_colours.mapped_colour= 127;
   end
   if ~isfield( eidors_colours, 'ref_level' );
      eidors_colours.ref_level= 0;
   end

   pp= eidors_colours;

function colours=set_mapped_colour(pp, scale, backgnd, elem_data)
   % need to generate a colourmap with pp.mapped_colour+1 elements
   % background pixel will be at entry #1. Thus for
   % mapped_colour= 3. CMAP = [backgnd,[-1 -.5  0 .5 1]
   %
   % Note: ensure patch uses 'direct' CDataMapping
   ncol= pp.mapped_colour;
   backgndidx= 1;
   [red,grn,blu] = blu_red_axis( pp, ...
          [-1,linspace(-1,1,2*ncol - 1)]', backgndidx );
   colormap([red,grn,blu]);
   colours = round( elem_data/ scale * (ncol-1))' + ncol + 1;
   colours(backgnd)= backgndidx;
