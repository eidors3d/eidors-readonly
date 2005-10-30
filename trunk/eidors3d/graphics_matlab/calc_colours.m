function colours= calc_colours(img, scale, do_colourbar)
% colours= calc_colours(img, scale, do_colourbar)
% Calculate a colour for each image element 
% Conductive (positive) areas are shown in red
% Non-Conductive (negative) areas are shown in blue
%
% img - an eidors image object
%    or
%     - a matrix of values;
%
% scale - colour value corresponding to maximum
%       - if not specified or scale==[] => autoscale
%
% do_colourbar ==1 => show a Matlab colorbar with
%        appropriate scaling
%
% Colour maps are controlled by the global variable
%   eidors_colours. The following settings
%
%   eidors_colours.greylev = .2; % Grey level to colour to
%   eidors_colours.sat_adj = .9; % max G,B when R=1
%   eidors_colours.backgnd= [.5,.5,.15]; % colour for non image regions
%   eidors_colours.mapped_colour= 0; % use colormap function
%         if mapped_colour is non-zero, it indicates the colourmap
%         size. Total colourmap is 2*mapped_colour

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_colours.m,v 1.13 2005-10-30 23:34:03 aadler Exp $  

% TODO: create a global eidors_colours object to control behaviour

pp=get_colours;
if nargin==0; return; end
 
if isfield(img,'type')
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

backgnd= isnan(elem_data);
elem_data(backgnd)= mean( elem_data(~backgnd));
e= length(elem_data);

% can't use | or || to support all stupid matlab versions > 6.0
if nargin <= 1
   scale =  max(abs(elem_data)) + eps;
elseif isempty(scale)
   scale =  max(abs(elem_data)) + eps;
end

if ~pp.mapped_colour
   [red,grn,blu] = blu_red_axis( pp, elem_data / scale, backgnd );
   colours= shiftdim( [red,grn,blu], -1);
else
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
end

% print colorbar if do_colourbar is specified
% can't use & or && to support all stupid matlab versions > 6.0
if nargin >=3 if do_colourbar==1
   hh= colorbar;
   p= get(hh,'Position');
   pm= p(2) + p(4)/2;
   set(hh,'Position', [p(1)+p(3), pm-p(4)*.6/2, p(3)*.6, p(4)*.6]);

   % set scaling
   OrdOfMag = 10^floor(log10(scale));
   scale_r  = OrdOfMag * floor( scale / OrdOfMag );
   ticks = 128.5 + 127.5*[-1,0,+1]*scale_r/scale;
   set(hh,'YTick', ticks');
   set(hh,'YTickLabel', [-scale_r, 0, scale_r]');
end; end


%scaled data must go from -1 to 1
function [red,grn,blu] = blu_red_axis( pp, scale_data, backgnd )
   D= sign(pp.greylev+eps); %force 0 to 1
   glev= abs(pp.greylev);
   F= 3*pp.sat_adj;

   red= D*F*abs(scale_data-1/F) - D + (D==-1);
   red= red.*(red>0).*(red<1) + (red>=1);
   red= red*(1-glev) + glev;

   grn= D*F*abs(scale_data    ) - D + (D==-1);
   grn= grn.*(grn>0).*(grn<1) + (grn>=1);
   grn= grn*(1-glev) + glev;

   blu= D*F*abs(scale_data+1/F) - D + (D==-1);
   blu= blu.*(blu>0).*(blu<1) + (blu>=1);
   blu= blu*(1-glev) + glev;

   red(backgnd) = pp.backgnd(1);
   grn(backgnd) = pp.backgnd(2);
   blu(backgnd) = pp.backgnd(3);


function pp=get_colours;
   global eidors_colours;

   if isempty( eidors_colours );
      eidors_colours.greylev = -.01;
      eidors_colours.sat_adj = .9;
      eidors_colours.backgnd= [.5,.5,.15];
      eidors_colours.mapped_colour= 0;
   end

   pp= eidors_colours;
