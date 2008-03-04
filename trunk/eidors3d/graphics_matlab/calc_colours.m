function colours= calc_colours(img, clim, do_colourbar, ref_lev)
% colours= calc_colours(img, clim, do_colourbar, ref_lev)
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
% Cs is 1xEx1 colourmap entries (if mapped_colour>0)
%       1xEx3 colourmap entries (if mapped_colour==0)
%
% Usage #2 (img is a MxN image matrix of reconstructed pixels):
%  
%   c_img = calc_colours( img);
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
%   'npoints' (DEFAULT 64) number of points accross the image
%   'clim'    (DEFAULT []) crop colour display of values above clim
% 
% PARAMETERS: clim
%    clim - colour limit. values more different from ref_level are cropped.
%         - if not specified or clim==[] => no limit
%    clim can be specified three ways (in decending priority order)
%       1. clim parameter to calc_colours('clim')
%       2. clim in img.calc_colours.clim
%       3. clim parameter to calc_colours()
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

% (C) 2005-2006 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_colours.m,v 1.46 2008-03-04 16:28:52 aadler Exp $  

if nargin==0;
% If no args - set defaults
    get_colours;
    return;
end

if isstr(img)
    % called as calc_colours('parameter' ... )
    test_exist_colours;
    if nargin==1;
       colours= get_field(img);
    else
       colours= set_field(img, clim);
    end
    return;

elseif isfield(img,'type')
   if strcmp( img.type, 'image' )
      try
         img_data= img.node_data;
      catch
         img_data= img.elem_data; %col vector
      end
   else
      error('calc_colours: input is not eidors image object');
   end
else
   img_data= img;
end

% Set default parameters
if nargin < 2; clim=[];                end
if nargin < 3; do_colourbar = 0;       end
if nargin < 4; ref_lev = 'use_global'; end


if isempty(img_data)
    colours = 'k'; %black
    return;
end

pp=get_colours;

m= size(img_data,1); n=size(img_data,2);

clim = calc_clim(img, clim);
% We can only plot the real part of data
% Vectorize img_data here, is get's reshaped later
[scl_data, ref_lev, max_scale] = ...
      scale_for_display( real(img_data(:)), ref_lev, clim );

backgnd= isnan(scl_data);
scl_data(backgnd)= mean( scl_data(~backgnd));

if pp.mapped_colour
   colours=set_mapped_colour(pp, backgnd, scl_data);
   if n==1;
      colours= reshape( colours, 1,[]);
   else
      colours= reshape( colours, m,n,[]);
   end
else
   [red,grn,blu] = blu_red_axis( pp, scl_data, backgnd );
   if n==1;
      colours= reshape( [red,grn,blu],1,[],3);
   else
      colours= reshape( [red,grn,blu],m,n,3);
   end
end

% print colorbar if do_colourbar is specified
if do_colourbar
   if ~pp.mapped_colour
       warning('Colorbar not available without mapped_colour option');
   else

   hh= colorbar; delete(hh); hh=colorbar;
   % make colourbar smaller and closer to axis
   p= get(hh,'Position');
   pm= p(2) + p(4)/2;
   set(hh,'Position', [p(1)+1.2*p(3), pm-p(4)*.6/2, p(3)*.6, p(4)*.6]);

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

   end
end


%scaled data must go from -1 to 1
function [red,grn,blu] = blu_red_axis( pp, scale_data, backgnd )
   % window data such that slope above w is 1/3 of that below
   % thus w is mapped to k st k/w = 3(1-k)/(1-w) -> k=3w/(1+2w)
   W= pp.window_range; K= 3*W/(1+2*W);
   scale_data= sign(scale_data) .* ( ...
     (  K/W*   abs(scale_data)   ) .* (abs(scale_data)<=W) + ...
     (K+K/W/3*(abs(scale_data)-W)) .* (abs(scale_data)> W) );

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

   red(backgnd) = pp.backgnd(1);
   grn(backgnd) = pp.backgnd(2);
   blu(backgnd) = pp.backgnd(3);


function test_exist_colours;
   global eidors_colours;

   if ~isfield( eidors_colours, 'greylev' );
      eidors_colours.greylev = -.001;
   end
   if ~isfield( eidors_colours, 'sat_adj' );
      eidors_colours.sat_adj = .9;
   end
   if ~isfield( eidors_colours, 'window_range' );
      eidors_colours.window_range = .5;
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
      eidors_colours.ref_level= 'auto';
   end
   if ~isfield( eidors_colours, 'npoints' );
      eidors_colours.npoints= 64;
   end
   if ~isfield( eidors_colours, 'clim' );
      eidors_colours.clim= [];
   end

function pp=get_colours;
   test_exist_colours;
   global eidors_colours;
   pp= eidors_colours;

function colours=set_mapped_colour(pp, backgnd, img_data)
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
   colours = fix( img_data * (ncol-1))' + ncol + 1;
   colours(backgnd)= backgndidx;

function value= get_field(param);
    global eidors_colours;
    value = getfield(eidors_colours, param);

function value= set_field(param, value);
    global eidors_colours;
    eidors_colours = setfield(eidors_colours, param, value);
    value= eidors_colours;

function clim = calc_clim(img, clim)
%       1. clim parameter to calc_colours('clim')
   if ~isempty(clim);
      return;
   end
%       2. clim in img.calc_colours.clim
   try
      clim= img.calc_colours.clim;
   catch;
%       3. clim parameter to calc_colours()
      global eidors_colours;
      clim = eidors_colours.clim;
   end
