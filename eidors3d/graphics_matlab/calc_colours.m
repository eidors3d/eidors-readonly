function colours= calc_colours(img, scale)
% colours= calc_colours(img)
% Calculate a colour for each image element 
% Conductive (positive) areas are shown in red
% Non-Conductive (negative) areas are shown in blue
%
% img - an eidors image object
%    or
%     - a matrix of values;
%
% scale - colour value corresponding to maximum
%       - default = autoscale
%
% Colour maps are controlled by the global variable
%   eidors_colours. The following settings
%
%   eidors_colours.greylev = .2; % Grey level to colour to
%   eidors_colours.sat_adj = .9; % max G,B when R=1
%   eidors_colours.backgnd= [.5,.5,.15]; % colour for non image regions
%   eidors_colours.mapped_colour= 0; % use colormap function
%         if mapped_colour is non-zero, it indicates the colourmap
%         size to use

% TODO: create a global eidors_colours object to control behaviour
% $Id: calc_colours.m,v 1.6 2005-10-27 12:02:00 aadler Exp $  

pp=get_colours;
 
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
[red,grn,blu] = blu_red_axis( pp, elem_data / scale, backgnd );
colours= ones(1, length(elem_data), 3);

glev= abs(pp.greylev);
   colours(1,:,:)= [red,grn,blu]*(1-glev) + glev;
if pp.greylev < 0
   colours = 1- colours(1,:,[3,2,1]);
end


function [red,grn,blu] = blu_red_axis( pp, scale_data, backgnd )
   F= 3*pp.sat_adj;
   red= F*abs(scale_data+1/F) -1;
   red= red.*(red>0).*(red<1) + (red>=1);
   red(backgnd) = pp.backgnd(1);

   grn= F*abs(scale_data    ) -1;
   grn= grn.*(grn>0).*(grn<1) + (grn>=1);
   grn(backgnd) = pp.backgnd(2);

   blu= F*abs(scale_data-1/F) -1;
   blu= blu.*(blu>0).*(blu<1) + (blu>=1);
   blu(backgnd) = pp.backgnd(3);


function pp=get_colours;
   global eidors_colours;

   if isempty( eidors_colours );
      eidors_colours.greylev = .2;
      eidors_colours.sat_adj = .9;
      eidors_colours.backgnd= [.5,.5,.15];
      eidors_colours.mapped_colour= 0;
   end

   pp= eidors_colours;
