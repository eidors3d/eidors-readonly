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
%   'mapped_colour' (DEFAULT 127) number of colourmap entries
%      using mapped_colour allows matlab to print vector graphics to eps
%      setting mapped_colour=0 means to use RGB colours
%   'npoints' (DEFAULT 64) number of points accross the image
%   'clim'    (DEFAULT []) crop colour display of values above clim
%           colour limit. values more different from ref_level are cropped.
%           if not specified or clim==[] => no limit
%   'cmap_type'  Specify special colours (Default 'blue_red')
%           if 'draeger' use the Draegerwerk/Amato colourmap
%           if 'jet' use the matlab jet colourmap
%   'cb_shrink_move' shrink or move the colorbar. See eidors_colourbar
%           help for details.
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

if nargin==0; error('calc_colours: expecting argument'); end
if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

if ischar(img)
    % called as calc_colours('parameter' ... )
    if nargin==1;
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
[scl_data, ref_lev, max_scale] = ...
      scale_for_display( real(img_data(:)), pp.ref_level, pp.clim );

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
   else
       eidors_colourbar(max_scale,ref_lev,pp.cb_shrink_move)
   end
end


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

   switch pp.cmap_type
     case {'blue_red',''}
      [red,grn,blu]= blue_red_colours(pp,scale_data);
     case 'draeger'
      [red,grn,blu]= draeger_colours(pp,scale_data);
     case 'jet'
      [red,grn,blu]= jet_colours(pp,scale_data);
     case 'jetair'
      % Make most of the range for air -ve
      scd = (scale_data + 0.5)/0.6;
      [red,grn,blu]= jet_colours(pp,scd);
     otherwise
      error(['specified cmap_type not understood:',pp.cmap_type]);
   end

   red(backgnd) = pp.backgnd(1);
   grn(backgnd) = pp.backgnd(2);
   blu(backgnd) = pp.backgnd(3);

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

function [red,grn,blu] = draeger_colours(pp,scale_data);
   grn=      (-scale_data>0.2) .* (-scale_data - 0.2)/0.8;
   red=grn + ( scale_data>0.2) .* ( scale_data - 0.2)/0.8*2/3;

   blu=      (-scale_data>0.1) .* (-scale_data - 0.1)/0.1;
   blu(-scale_data>0.2) = 1;
   blu=blu + ( scale_data>0.2) .* ( scale_data - 0.2)/0.8*2/3;

% Sometimes this is just slightly > 1
   red(red>1) = 1;
   grn(grn>1) = 1;
   blu(blu>1) = 1;

function [red,grn,blu] = jet_colours(pp,scale_data);
   grn = 2*(3/4 - abs(scale_data));
   grn(grn>1) = 1; grn(grn<0) = 0;

   red = 1.5 - 2*abs(scale_data + 0.5);
   red(red>1) = 1; red(red<0) = 0;

   blu = 1.5 - 2*abs(scale_data - 0.5);
   blu(blu>1) = 1; blu(blu<0) = 0;

function pp=get_colours( img );
   global eidors_colours;
   pp= eidors_colours;

% override global if calc.colours specified
   try
% DAMN Matlab should have syntax for this loop
      fds= fieldnames(img(1).calc_colours);%assume all are the same!
      for fdn= fds(:)';
         fdn= fdn{1};
         pp.( fdn ) = img(1).calc_colours.(fdn);
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
   img = eidors_obj('image','test'); 

   img.calc_colours.mapped_colour = 127;
   img.calc_colours.ref_level = 'auto';
   img.calc_colours.sat_adj = 0.9;
   img.calc_colours.window_range = 0.9;
   img.calc_colours.greylev = -0.01;
 
   img.elem_data = [-2;0;0;0;1;3];
   unit_test_cmp('cc01', calc_colours(img), [44; 128; 128; 128; 170; 254]);

   imgk(1) = img; imgk(2) = img;
   unit_test_cmp('cc01', calc_colours(imgk), [44; 128; 128; 128; 170; 254]*[1,1]);

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
   unit_test_cmp('cc05',cc, [2,2,2;4,4,4;2,4,4;0,1,4;0,0,0;4,1,0;4,4,2;4,4,4]/4,1e-4);

% TESTS TO WRITE
%   'greylev'    (DEFAULT -.01)
%   'sat_adj'    (DEFAULT .9)
%   'window_range' (DEFAULT .9)
%   'backgnd' ( DEFAULT [.5,.5,.15] )
%   'ref_level' (DEFAULT 'auto')
%   'mapped_colour' (DEFAULT 127)
%   'npoints' (DEFAULT 64)
%   'clim'    (DEFAULT [])
%   'cmap_type'  (Default blue-red)

