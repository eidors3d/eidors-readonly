function [colours,scl_data]= calc_colours(img, set_value, do_colourbar)
% [colours,scl_data]= calc_colours(img, set_value, do_colourbar)
% Calculate a colour for each image element 
%
% Conductive (positive) areas are shown in red
% Non-Conductive (negative) areas are shown in blue
%
% PARAMETERS: img
%     - 1) an EIDORS image object, or a Ex1 vector
%     - 2) a 2D image matrix
%
% Usage #1 (img is a Ex1 vector of element conductivities)
%
%   Cs = calc_colours( img);
%   patch(Xs,Ys,Zs,Cs);
%
% Cs is Ex1x1 colourmap entries (if mapped_colour>0)
%       Ex1x3 colourmap entries (if mapped_colour==0)
%
% Usage #2 (rimg is a MxN image matrix of reconstructed pixels):
%           img is an image structure with the image properties
%  
%   c_img = calc_colours( rimg, img);
%   image( c_img );
%
% c_img is MxN colourmap entries (if mapped_colour>0)
%          MxNx3 colourmap entries (if mapped_colour==0)
%
% Usage #3 (img is string parameter value)
%  
%   value = calc_colours( 'param' );
%   calc_colours( 'param', value );
%    eg. calc_colours('mapped_colour',127)
%         Use this to allow printing vector eps files to get around
%         a matlab bug with warning 'RGB color data not ...'
%
%   calc_colours( 'component', 'real' ); (DEFAULT real)
%     the other value is 'imag' to show the imaginary component
%
%   The following parameters are accepted
%
%   'greylev'    (DEFAULT -.01): the colour of the ref_level.
%      Negative values indicate white background
%      For almost white, greylev=-.01; Black=> greylev=.01
%   'sat_adj'    (DEFAULT .9): max G,B when R=1
%   'window_range' (DEFAULT .9); window colour range
%      Colour slope outside range is 1/3 of centre slope
%   'backgnd' ( DEFAULT [.5,.5,.15] ): image border colour 
%   'ref_level' (DEFAULT 'auto') conductivity of centre of
%      colour mapping. 'auto' tries to estimate a good level.
%      For complex image data, ref_level should also be complex.
%   'mapped_colour' (DEFAULT 127) number of colourmap entries
%      using mapped_colour allows matlab to print vector graphics to eps
%      setting mapped_colour=0 means to use RGB colours
%   'npoints' (DEFAULT 64) number of points accross the image
%   'transparency_thresh' fraction of maximum value at which colours
%             are rendered transparent for 
%   'clim'    (DEFAULT []) crop colour display of values above clim
%           colour limit. values more different from ref_level are cropped.
%           if not specified or clim==[] => no limit
%   'cmap_type'  Specify special colours (Default 'blue_red')
%           'blue_red':          default Blue/Red eidors colourmpa
%           'jet':               matlab jet colourmap
%           'jetair':            scaled jet colours
%           'blue_yellow':       Blue/Yellow colours
%           'greyscale':         Greyscale colours (Lungs white)
%           'greyscale-inverse': Greyscale colours (Lungs black)
%           'copper':            Copper colours
%           'blue_white_red':    Blue/White/Red colours
%           'black_red':         Black/Red Colours
%           'blue_black_red':    Blue/Black/Red colours
%           'polar_colours':     "Polar" blue/white/red colours
%           'draeger', 'draeger-difference':
%                                Draegerwerk colourmap (difference)
%           'draeger-2009':      Draegerwerk colourmap (2009)
%           'draeger-tidal':     Draegerwerk colourmap (tidal images)
%           'swisstom'           Swisstom colourmap
%           'timpel'             Timpel colourmap
%           'flame'              White/Red/Yellow/Blue colours
%           'ice'                White/Purple/Black/Orange/Yellow colours
%           matrix [Nx3]         Use values ([0-1]) as REB colourmap
%   'cb_shrink_move' shrink or move the colorbar. See eidors_colourbar
%           help for details.
%   'image_field', 'image_field_idx', 'image_field_val' 
%           image elements which match _idx are set to _val.
%   'defaults' set to default colours in list above
%
%   'colourmap' Return the current EIDORS colormap. 
%           Use as colormap(calc_colours('colourmap'))
%
% PARAMETERS CAN BE SPECIFIED IN TWO WAYS
%   1. as an image parameter (ie clim in img.calc_colours.clim)
%   2. a second parameter to ( calc_colours(data, param2 )
%          where param2.calc_colours.clim= ... etc
%   3. parameter to calc_colours('clim')
%
% Parameters specified as (1) will override (2)
%
% PARAMETERS: do_colourbar
%    - show a Matlab colorbar with appropriate scaling
%
%  usage: c_img= calc_colours( img, clim );
%         image( c_img );
%         calc_colours( img, clim, 1); %now do colorbar 
%
% PARAMETERS: ref_lev
%     - if specified, override the global ref_level parameter
%

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id$  

% BUGS TO FIX:
%  mapped_colour should set the colourmap length for all

if nargin==0; error('calc_colours: expecting argument'); end
if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

if ischar(img)
    % called as calc_colours('parameter' ... )
    if nargin==1;
       if strcmp(img,'defaults') % set defaults and return
          set_colours_defaults; return
       end
       colours= get_field(img);
    else
       colours= set_field(img, set_value);
    end
    return;

elseif isfield(img,'type')
   img_data= get_img_data( img );
   pp=get_colours(img(1)) ;
else
   img_data= img;

   if nargin==1
      pp=get_colours( [] ); 
   else
      pp=get_colours(set_value); 
   end
end

% Set default parameters
if nargin < 3; do_colourbar = 0;       end
ref_lev = 'use_global';


if isempty(img_data)
    colours = 'k'; %black
    return;
end

m= size(img_data,1); n=size(img_data,2);

% We can only plot the real part of data
% Vectorize img_data here, it gets reshaped later
[scl_data, ref_lev, max_scale] = scale_for_display( img_data(:), pp);

backgnd= isnan(scl_data);
scl_data(backgnd)= mean( scl_data(~backgnd));

if pp.mapped_colour
   colours=set_mapped_colour(pp, backgnd, scl_data);
   colours= reshape( colours, m,n,[]);
else
   [red,grn,blu] = blu_red_axis( pp, scl_data, backgnd );
   colours= reshape( [red,grn,blu],m,n,3);
end

% print colorbar if do_colourbar is specified
if do_colourbar
   if ~pp.mapped_colour
       warning('Colorbar not available without mapped_colour option');
   elseif do_colourbar < 0
        eidors_colourbar(-do_colourbar,max_scale,pp.cb_shrink_move,1)
   else
       eidors_colourbar(max_scale,ref_lev,pp.cb_shrink_move)
   end
end


function set_colours_defaults;
   clear -global eidors_colours;

   calc_colours('greylev',-.001);          % background colour = white
   calc_colours('sat_adj',.9);             % saturation of red and blue
   calc_colours('window_range', .7);       % windowing of colours
   calc_colours('backgnd',[.35,.50,.50]);    % background colour
   calc_colours('mapped_colour',127);      % use 127*2+1 colourmap entries
   calc_colours('ref_level','auto');       % auto set background colour
   calc_colours('npoints',64);             % 64 raster points
   calc_colours('clim',[]);                % no colour cropping
   calc_colours('cb_shrink_move',[1,1,0]); % Don't shrink or move colorbar
   calc_colours('transparency_thresh',0.25); % transparent at .25 of max
   eidors_msg('Setting Default Colours',2);
   calc_colours('cmap_type','blue_red');   % default eidors colours


%scaled data must go from -1 to 1
function [red,grn,blu] = blu_red_axis( pp, scale_data, backgnd )
   if ~isfield(pp,'cmap_type')
      pp.cmap_type = 'blue_red';
   end; 

   % window data such that slope above w is 1/3 of that below
   % thus w is mapped to k st k/w = 3(1-k)/(1-w) -> k=3w/(1+2w)
   W= pp.window_range; K= 3*W/(1+2*W);
   scale_data= sign(scale_data) .* ( ...
     (  K/W*   abs(scale_data)   ) .* (abs(scale_data)<=W) + ...
     (K+K/W/3*(abs(scale_data)-W)) .* (abs(scale_data)> W) );

   if isnumeric(pp.cmap_type)
      user_colourmap = pp.cmap_type;
      pp.cmap_type= 'user_colourmap';
   end

   switch pp.cmap_type
     case 'user_colourmap'
      red = user_colourmap(:,1);
      grn = user_colourmap(:,2);
      blu = user_colourmap(:,3);
     case 'blue_red'
      [red,grn,blu]= blue_red_colours(pp,scale_data);
     case {'draeger', 'draeger-difference'}
      [red,grn,blu]= draeger_colours(pp,scale_data, 'difference2014');
     case 'draeger-2009',
      [red,grn,blu]= draeger_colours(pp,scale_data, '2009');
     case 'draeger-tidal',
      [red,grn,blu]= draeger_colours(pp,scale_data, 'tidal2014');
     case 'swisstom'
      [red,grn,blu]= swisstom_colours(pp,scale_data);
     case 'timpel'
      [red,grn,blu]= timpel_colours(pp,scale_data);
     case 'jet'
      [red,grn,blu]= jet_colours(pp,scale_data);
     case 'jetair'
      % Make most of the range for air -ve
      scd = (scale_data + 0.5)/0.6;
      [red,grn,blu]= jet_colours(pp,scd);
     case 'blue_yellow'
         [red,grn,blu] = blue_yellow_colours(pp,scale_data);
     case 'greyscale'          % Lungs are white
         [red,grn,blu] = greyscale_colours(pp,scale_data);
     case 'greyscale-inverse'  % Lungs are black
         [red,grn,blu] = greyscale_colours(pp,-scale_data);
     case 'copper'
         [red,grn,blu] = copper_colours(pp,scale_data);
     case 'blue_white_red'
         [red,grn,blu]= blue_white_red_colours(pp,scale_data);
     case 'black_red'
         [red,grn,blu] = black_red_colours(pp,scale_data);
     case 'blue_black_red'
         [red,grn,blu] = blue_black_red_colours(pp,scale_data);
     case 'polar_colours'
         [red,grn,blu] = polar_blue_white_red_colours(pp,scale_data);
     case 'flame'
         [red,grn,blu] = flame_colours(pp,scale_data);
     case 'ice'
         [red,grn,blu] = ice_colours(pp,scale_data);

     otherwise
      error(['specified cmap_type not understood:',pp.cmap_type]);
   end
% Sometimes this is just slightly > 1
   red(red>1) = 1;
   grn(grn>1) = 1;
   blu(blu>1) = 1;


   red(backgnd) = pp.backgnd(1);
   grn(backgnd) = pp.backgnd(2);
   blu(backgnd) = pp.backgnd(3);

%FIXME - come up with a general way to handle blocks in image field
if isfield(pp,'image_field')
    red(2)=pp.image_field(1);
    grn(2)=pp.image_field(2);
    blu(2)=pp.image_field(3);
end

   
function [red,grn,blu]= blue_red_colours(pp,scale_data);
if 0
   D= sign(pp.greylev+eps); %force 0 to 1
   glev= abs(pp.greylev);
   F= 3*pp.sat_adj;

   red= D*F*abs(scale_data+D/F) - D + (D==-1);
   red= red.*(red>0).*(red<1) + (red>=1);
   red= red*(1-glev) + glev;
end
   ofs= (pp.greylev >= 0);   % 1 if greylev>=0
   glev= abs(pp.greylev);

   D= (2*ofs - 1);
   ofs= ofs - 2*(ofs==0);
   F= 3*pp.sat_adj;
   DF= D*F; D_F= D/F;

   red= DF*abs(scale_data+D_F) - ofs;
   red= red.*(red>0).*(red<1) + (red>=1);

   grn= DF*abs(scale_data    ) - ofs;
   grn= grn.*(grn>0).*(grn<1) + (grn>=1);

   blu= DF*abs(scale_data-D_F) - ofs;
   blu= blu.*(blu>0).*(blu<1) + (blu>=1);

   if pp.greylev >=0 % Black background
      red= red*(1-glev) + glev;
      grn= grn*(1-glev) + glev;
      blu= blu*(1-glev) + glev;
   else
      red= red*(1-glev);
      grn= grn*(1-glev);
      blu= blu*(1-glev);
   end

function [red,grn,blu]= swisstom_colours(pp,scale_data);
   backgnd =  hex2dec(['59';'59';'59']);
   alphapts = [-255,32,160,255]/255;
   alphaval = [0, 0, 1, 1];
   alp = interp1(alphapts, alphaval, -scale_data, 'linear');

   hexcolours = ['3E4A50';
                 '2F4D58';
                 '357892';
                 '2AB1DF';
                 '61D1FB';
                 'CFF0FE'];
   redpts = hex2dec(hexcolours(:,1:2));
   grnpts = hex2dec(hexcolours(:,3:4));
   blupts = hex2dec(hexcolours(:,5:6));
   interppts = linspace(0,1, size(hexcolours,1));
% add for +ve changes
   redpts = redpts([1,1:end]);
   grnpts = grnpts([1,1:end]);
   blupts = blupts([1,1:end]);
   interppts = [-1,interppts];

   red = interp1(interppts, redpts,-scale_data, 'linear');
   grn = interp1(interppts, grnpts,-scale_data, 'linear');
   blu = interp1(interppts, blupts,-scale_data, 'linear');

   red = red.*alp + backgnd(1)*(1-alp);
   grn = grn.*alp + backgnd(2)*(1-alp);
   blu = blu.*alp + backgnd(3)*(1-alp);
          
   red = red/255;
   grn = grn/255;
   blu = blu/255;

function [red,grn,blu]= timpel_colours(pp,scale_data)
   sd  = -scale_data*100;
   P00 = (sd<  9);
   p09 = (sd>= 9) .* (sd<10) .* (sd- 9) / (10- 9); P09 = p09 >0;
   p10 = (sd>=10) .* (sd<17) .* (sd-10) / (17-10); P10 = p10 >0;
   p17 = (sd>=17) .* (sd<30) .* (sd-17) / (30-17); P17 = p17 >0;
   p30 = (sd>=30) .* (sd<60) .* (sd-30) / (60-30); P30 = p30 >0;
   p60 = (sd>=60) .* (sd<85) .* (sd-60) / (85-60); P60 = p60 >0;
   p85 = (sd>=85)            .* (sd-85) /(100-85); P85 = p85 >0;


   red =  150*P00 + 198*p85;
   grn =  150*P00 +       88 *P09 + ...
           88*P10 +  (99- 88)*p10 + ...
           99*P17 + (121- 99)*p17 + ...
          121*P30 + (190-121)*p30 + ...
          190*P60 + (243-190)*p60 + ...
          243*P85 + (251-243)*p85;
   blu =  150*P00 +      159 *P09 + ...
          180*P10 + (180-159)*p10 + ...
          255*P17 + ...
          255*P30 + ...
          255*P60 + ...
          255*P85 ;
          
   red = red/255;
   grn = grn/255;
   blu = blu/255;

% From Vinicius Torsani
%These are the color codes from range of 0 to 1:
% 
%         R       G       B
% 0 to 0,09       150     150     150
% to 0,1  0       88      159
% to 0,17         0       99      180
% to 0,3  0       121     255
% to 0,6  0       190     255
% to 0,85         0       243     255
% to 1    198     251     255




function [red,grn,blu] = draeger_colours(pp,scale_data, version);
   switch version
     case '2009'
       grn=      (-scale_data>0.2) .* (-scale_data - 0.2)/0.8;
       red=grn + ( scale_data>0.2) .* ( scale_data - 0.2)/0.8*2/3;
       
       blu=      (-scale_data>0.1) .* (-scale_data - 0.1)/0.1;
       blu(-scale_data>0.2) = 1;
       blu=blu + ( scale_data>0.2) .* ( scale_data - 0.2)/0.8*2/3;
     case 'tidal2014'
       sd  = -scale_data*100;
       p20 = (sd> 20) .* (sd<= 100) .* (sd - 20) / (100-20);
       P20 = (p20>0);
       p10 = (sd> 10) .* (sd<=  20) .* (sd - 10) / ( 20-10);
       n20 = (sd<-20) .* (sd>=-100) .* (sd + 20) /-(100-20);
       red = 255*p20 + 170*n20;
       grn = 255*p20;
       blu = 255*P20 + 255*p10 + 170*n20;

       red = red/255;
       grn = grn/255;
       blu = blu/255;

     case 'difference2014';
       sd  = -scale_data*100;
       p40 = (sd> 40) .* (sd<= 100) .* (sd - 40) / (100-40);
       P40 = (p40>0);
       p10 = (sd> 10) .* (sd<=  40) .* (sd - 10) / ( 40-10);
       n10 = (sd<-10) .* (sd>=- 40) .* (sd + 10) /-( 40-10);
       n40 = (sd<-40) .* (sd>=-100) .* (sd + 40) /-(100-40);
       N40 = (n40>0);
       red =  132*p40  + 255*n10  + 255*N40;
       grn =  153*p10  + 153*P40 + (240-153)*P40 + ...
              153*n10  + 153*N40 + (240-153)*n40;
       blu =  153*n40  + 255*p10 +  255*P40;

       red = red/255;
       grn = grn/255;
       blu = blu/255;
     otherwise; error('version "%s" not recognized', version);
   end


function [red,grn,blu] = jet_colours(pp,scale_data);
   grn = 2*(3/4 - abs(scale_data));
   grn(grn>1) = 1; grn(grn<0) = 0;

   red = 1.5 - 2*abs(scale_data + 0.5);
   red(red>1) = 1; red(red<0) = 0;

   blu = 1.5 - 2*abs(scale_data - 0.5);
   blu(blu>1) = 1; blu(blu<0) = 0;

% TODO: only works with mapped_colours, fix to use scale_data
function [red,grn,blu] = blue_yellow_colours(pp,scale_data);
   % Inpired by 
   % FIREICE LightCyan-Cyan-Blue-Black-Red-Yellow-LightYellow Colormap
   % by: Joseph Kirk %  jdkirk630@gmail.com
   clrs = flipud([0.75 1 1; 0 1 1; 0 0 1;...
       0 0 0; 1 0 0; 1 1 0; 1 1 0.75]);

   y = linspace(-1,1,size(clrs,1));
   cc = interp2(1:3,y,clrs,1:3,scale_data);
   red = cc(:,1);
   grn = cc(:,2);
   blu = cc(:,3);


% TODO: only works with mapped_colours, fix to use scale_data
function [red,grn,blu]= greyscale_colours(pp,scale_data);
   cc= 0.5 - scale_data/2;
   red = cc;
   grn = cc;
   blu = cc;

   

% TODO: only works with mapped_colours, fix to use scale_data
function [red,grn,blu]= copper_colours(pp,scale_data);
% Constants from the copper.m function
   cc= 0.5 + scale_data/2;
   red = min(cc * 1.2500,1);
   grn = cc * 0.7812;
   blu = cc * 0.4975;

function [red,grn,blu]= blue_white_red_colours(pp,scale_data);
   % bluewhitered Colormap
   % Based on FIREICE from Joseph Kirk jdkirk630@gmail.com
   clrs = flipud([0 0 1; 1 1 1; 1 0 0]);

   y = linspace(-1,1,size(clrs,1));
   cc = interp2(1:3,y,clrs,1:3,scale_data);
   red = cc(:,1);
   grn = cc(:,2);
   blu = cc(:,3);



function [red,grn,blu] = black_red_colours(pp,scale_data);
   red = 0.5 + scale_data/2;
   blu = 0*scale_data;
   grn = blu;

function [red,grn,blu] = flame_colours(pp,scale_data);
   sd = 0.5 + scale_data/2; % 0..1
   red = sd*0+1;
   blu = (1 - 3*sd);
   grn = (1 - 3*sd);% + (sd>t).*(sd-t)/t;
   grn(grn<0) = -2*grn(grn<0);
   t=0.8;
   grn(sd>t) = 1-(sd(sd>t)-t)/(1-t);
   red(sd>t) = 1-(sd(sd>t)-t)/(1-t);
   blu(sd>t) = (sd(sd>t)-t)/(1-t);
   [red,grn,blu] = saturate(red,grn,blu);

function [red,grn,blu] = ice_colours(pp,scale_data);
   sd = 0.5 + scale_data/2; % 0..1
   blu = (1 - 2*sd);
   grn = (1 - 3*sd);
   red = (1 - 2*sd);
   red(red<0) = -2*red(red<0);
   t=0.5;
%   red(sd>t) = 1-(sd(sd>t)-t)/(1-t);
%   blu(sd>t) = 1-(sd(sd>t)-t)/(1-t);
   grn(sd>t) = (2*sd(sd>t)-1);
   [red,grn,blu] = saturate(red,grn,blu);

function [red,grn,blu] = saturate(red,grn,blu);
   red(red>1)=1;
   grn(grn>1)=1;
   blu(blu>1)=1;
   red(red<0)=0;
   grn(grn<0)=0;
   blu(blu<0)=0;

% TODO: only works with mapped_colours, fix to use scale_data
function [red,grn,blu] = blue_black_red_colours(pp,scale_data);
   red =-min(0,scale_data);
   blu = max(0,scale_data);
   grn = 0*scale_data;

% TODO: only works with mapped_colours, fix to use scale_data
% Ideas for this colourmap, and a few lines of code are from
% (C) 2011, Francois Beauducel, IPGP in polarmap.m
function [red,grn,blu] = polar_blue_white_red_colours(pp,scale_data);
%  0   0   1
%  .5  .5  1
%  1   1   1
%  1   .5  .5
%  1   0   0
   red = min(1, 1 + scale_data);
   grn = 1 - abs(scale_data);
   blu = min(1, 1 - scale_data);


function pp=get_colours( img );
   global eidors_colours;
   pp= eidors_colours;

   pp.component = 'real'; % Don't get from globals

% override global if calc.colours specified
   try
% DAMN Matlab should have syntax for this loop
      fds= fieldnames(img(1).calc_colours);%assume all are the same!
      for fdn= fds(:)';
         fdn= fdn{1};
         pp.( fdn ) = img(1).calc_colours.(fdn);
      end
   end
   if isfield(img, 'eidors_colourbar')
      pp.cb_shrink_move = [1,1,0];
      try
         pp.cb_shrink_move = img.eidors_colourbar.cb_shrink_move;
      end
   end

function BGI= backgndidx;
  BGI= 1;

function colours=set_mapped_colour(pp, backgnd, img_data)
   % need to generate a colourmap with pp.mapped_colour+1 elements
   % background pixel will be at entry #1. Thus for
   % mapped_colour= 3. CMAP = [backgnd,[-1 -.5  0 .5 1]
   %
   % Note: ensure patch uses 'direct' CDataMapping
   ncol= pp.mapped_colour;
   [red,grn,blu] = blu_red_axis( pp, ...
          [-1,linspace(-1,1,2*ncol - 1)]', backgndidx );
   colormap([red,grn,blu]);

   %FIXME - come up with a general way to handle blocks in image field
   if isfield(pp,'image_field')
      warning(['calc_colours: image_field is an experimental feature ' ...
                              'and will be re-implemented']);
      infmask = isinf(img_data);
      mindata = min(img_data(~infmask));
      colours = fix( (img_data-mindata) * (ncol-1))' + 3;
      colours(backgnd)= backgndidx;
      colours(pp.image_field_idx)=2;
      return
   end

% This is the line I wan't to write. However, matlab seems to waste
%  lots of memory to do it. Instead we break it into pieces
%  colours = fix( img_data * (ncol-1))' + ncol + 1;
   colours = img_data'; colours = colours * (ncol-1);
   colours = fix ( colours ); colours = colours + ncol + 1;
   colours(backgnd)= backgndidx;

function value= get_field(param);
   global eidors_colours;
   if strcmp(param,'colourmap')
      ncol= eidors_colours.mapped_colour;
      [red,grn,blu] = blu_red_axis( eidors_colours, ...
             [-1,linspace(-1,1,2*ncol - 1)]', backgndidx );
      value= [red,grn,blu];
   else
      value = getfield(eidors_colours, param);
   end

function value= set_field(param, value);
   global eidors_colours;
   eidors_colours = setfield(eidors_colours, param, value);
   value= eidors_colours;

% TESTS:
function do_unit_test
   calc_colours('defaults');

   img = eidors_obj('image','test'); 

   img.calc_colours.mapped_colour = 127;
   img.calc_colours.ref_level = 'auto';
   img.calc_colours.sat_adj = 0.9;
   img.calc_colours.window_range = 0.9;
   img.calc_colours.greylev = -0.01;
 
   img.elem_data = [-2;0;0;0;1;3];
   unit_test_cmp('cc01a', calc_colours(img), [44; 128; 128; 128; 170; 254]);

   imgk(1) = img; imgk(2) = img;
   unit_test_cmp('cc01b', calc_colours(imgk), [44; 128; 128; 128; 170; 254]*[1,1]);
   img.calc_colours.ref_level = 0;
   unit_test_cmp('cc01c', calc_colours(img), [44; 128; 128; 128; 170; 254]);

   img.calc_colours.ref_level = 1;
   unit_test_cmp('cc02', calc_colours(img), [ 2;  86;  86;  86; 128; 212]);

   img.calc_colours.greylev = -.1;
   img.calc_colours.mapped_colour = 0;
   img.elem_data = [-2;1;3];
   cm= reshape([ 0, 0.9, 0.9; 0, 0.9, 0.0643; 0.27, 0.9, 0]',[3,1,3]);
   unit_test_cmp('cc03', calc_colours(img), cm,1e-3)

   img.calc_colours.greylev =  .1;
   cm= reshape([ 0.73, 1.0, 1.0; 0.1, 0.1, 0.1; 1.0, 0.9357, 0.1],[3,1,3]);
   unit_test_cmp('cc04', calc_colours(img), cm,1e-3)

   calc_colours('greylev',0);
   calc_colours('sat_adj',1);
   calc_colours('mapped_colour',4);
   calc_colours('backgnd',[0.5,0.5,0.5]);
   cc= calc_colours('colourmap');
   unit_test_cmp('cc05',cc, ...
      [2,2,2;4,4,4;2,4,4;0,1,4;0,0,0;4,1,0;4,4,2;4,4,4]/4,1e-10);

   calc_colours('greylev',-.01);
   cc= calc_colours('colourmap');
   unit_test_cmp('cc06',cc, [200,200,200;0,0,0;0,0,198; ...
        0,297,396;396,396,396;396,297,0;198,0,0;0,0,0]/400,1e-10);

   calc_colours('greylev',0);
   calc_colours('sat_adj',.5);
   cc= calc_colours('colourmap');
   unit_test_cmp('sa01',cc(2,:), [0,5,10]/10,1e-10);

   calc_colours('sat_adj',.8);
   cc= calc_colours('colourmap');
   unit_test_cmp('sa02',cc(2,:), [4,10,10]/10,1e-10);

   calc_colours('sat_adj',1 );
   cc= calc_colours('colourmap');
   unit_test_cmp('sa03',cc(2,:), [10,10,10]/10,1e-10);

%   'window_range' (DEFAULT .9)
   unit_test_cmp('bg01',cc(1,:), [50,50,50]/100, 1e-10);
   calc_colours('defaults');
   cc= calc_colours('colourmap');
%  unit_test_cmp('bg01',cc(1,:), [50,50,15]/100, 1e-10); %EIDORS 3.6
%  unit_test_cmp('bg02',cc(1,:), [20,40,40]/100, 1e-10); %EIDORS 3.7
%  unit_test_cmp('bg02',cc(1,:), [35,45,45]/100, 1e-10); %EIDORS 3.8
   unit_test_cmp('bg02',cc(1,:), [35,50,50]/100, 1e-10); %EIDORS 3.9

   unit_test_cmp('mc01',size(cc), [127*2,3]);
   calc_colours('mapped_colour',4);
   cc= calc_colours('colourmap');
   unit_test_cmp('mc02',size(cc), [4*2,3]);
%   'ref_level' (DEFAULT 'auto')
%   'clim'    (DEFAULT [])
%   'cmap_type'  (Default blue-red)

   test.blue_red = [0,0,2997; 9990,9990,9990; 2997,0,0]/1e4;
%  test.draeger  = [10000,10000,10000;0,0,0;6667,0,6667]/1e4; % OLD 2009!
   test.draeger  = [5176,9412,10000;0,0,0;10000,9412,6000]/1e4; % 2014
   test.timpel   = [7765,9843,10000;5882,5882,5882;5882,5882,5882]/1e4;
   test.jet      = [5000,0,0;5000,10000,5000;0,0,5000]/1e4;
   test.jetair   = [8333,0,0; 0,0,8333;0,0,0]/1e4;
   test.blue_yellow = [10000,10000,7500; 0,0,0;7500,10000,10000]/1e4;
   test.greyscale = flipud([0,0,0;5000,5000,5000;10000,10000,10000])/1e4;
   test.copper = [0,0,0;6250,3906,2487;10000,7812,4975]/1e4;
   test.blue_white_red = [10000,0,0;10000,10000,10000;0,0,10000]/1e4;
   test.black_red = [0,0,0;5000,0,0;10000,0,0]/1e4;
   test.blue_black_red = [10000,0,0;0,0,0;0,0,10000]/1e4;
   test.polar_colours = [0,0,10000;10000,10000,10000;10000,0,0]/1e4;
   test.flame = [10000,10000,10000;10000,10000,0;0,0,10000]/1e4;
   test.ice = [10000,10000,10000;0,0,0;10000,10000,0]/1e4;

   calc_colours('defaults');
   calc_colours('mapped_colour',4);
   for ct = fieldnames(test)'
      ct = ct{1};
      tname = sprintf('ctm4(%s)',ct);
      calc_colours('cmap_type', ct);
      cc= calc_colours('colourmap');
      unit_test_cmp(tname,cc([2,5,8],:), test.( ct ), 1e-4);
   end

   calc_colours('defaults');
   calc_colours('mapped_colour',0);
   for ct = fieldnames(test)'
      ct = ct{1};
      tname = sprintf('ctm0(%s)',ct);
      calc_colours('cmap_type', ct);

      cc= calc_colours(linspace(-1,1,7)); cc= squeeze(cc);
      unit_test_cmp(tname,cc([1,4,7],:), test.( ct ), 1e-4);
   end

   calc_colours('defaults');
