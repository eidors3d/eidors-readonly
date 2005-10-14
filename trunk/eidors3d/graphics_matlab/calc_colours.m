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

% TODO: create a global eidors_colours object to control behaviour
% $Id: calc_colours.m,v 1.3 2005-10-14 18:27:26 aadler Exp $  

greylev=.2;
sat_adj=.9 ;
back_gnd_clr= [.5,.5,.15];
 
if isfield(img,'type')
   if strcmp( img.type, 'image' )
      elem_data= img.elem_data(:); %col vector
   else
      error('calc_colours: input is not eidors image object');
   end
else
   elem_data= img(:);
end

backgnd= isnan(elem_data);
elem_data(backgnd)= mean( elem_data(~backgnd));
e= length(elem_data);

% can't use | or || to support all stupid matlab versions > 6.0
if nargin <= 1
   scale =  max(abs(elem_data));
elseif isempty(scale)
   scale =  max(abs(elem_data));
end
scale_ed = elem_data / scale;

F= 3*sat_adj;
red= F*abs(scale_ed+1/F) -1;
red= red.*(red>0).*(red<1) + (red>=1);
grn(backgnd) = back_gnd_clr(1);

grn= F*abs(scale_ed    ) -1;
grn= grn.*(grn>0).*(grn<1) + (grn>=1);
grn(backgnd) = back_gnd_clr(2);

blu= F*abs(scale_ed-1/F) -1;
blu= blu.*(blu>0).*(blu<1) + (blu>=1);
grn(backgnd) = back_gnd_clr(3);

colours= ones(1, length(elem_data), 3);
colours(1,:,:)= [red,grn,blu]*(1-greylev)+ greylev;
