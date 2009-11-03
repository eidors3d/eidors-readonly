function params = eval_GREIT_fig_merit(imgs, xyzr_pt)
% EVAL_GREIT_FIG_MERIT: calculate GREIT figures of merit for images
% USAGE:
% params = eval_GREIT_fig_merit(imgs, xyzr_pt)
%  params(1,:) = Image Amplitude
%  params(2,:) = Position Error => + toward centre, - toward edge
%  params(3,:) = Resolution
%  params(4,:) = Shape Deformation
%  params(5,:) = Ringing
%
%  imgs:    a sequence of eidors images of single point targets
%  xyzr_pt: [x;y;z;radius] of each images point

% (C) 2009 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% TODO: make this more generic for different grids
imgs.calc_colours.npoints = 32;  % Eval on 32x32 grid
imgs = calc_slices(imgs);
imgs(isnan(imgs)) = 0;

N_imgs = size(imgs,3);
for i= 1:N_imgs
   [xmean,ymean,equiv_circ,map,qmi,img] = calc_cofg(imgs(:,:,i));
   params(1,i) = calc_amplitude( img );
   params(2,i) = calc_posn_error( qmi, xmean, ymean, xyzr_pt(1:2,i) );
   params(3,i) = calc_resolution( qmi, map );
   params(4,i) = calc_shape_deform( qmi, equiv_circ );
   params(5,i) = calc_ringing( img, qmi );
end

% TODO: Fix this when we start to care about units
params(1,:) = params(1,:)/mean(params(1,1:10));

function ampl = calc_amplitude(img)
   ampl = sum(img(:));

function pe   = calc_posn_error(qmi, xmean, ymean, xy)
   pe = sqrt(sum(xy.^2)) - sqrt( xmean^2 + ymean^2);

function res  = calc_resolution(qmi, map)
   res = sqrt( sum(qmi(:)) / sum(map(:)));

function sd  = calc_shape_deform(qmi, equiv_circ)
   not_circ= qmi & ~equiv_circ;
   sd = sum(not_circ(:))/sum(qmi(:));

function rr = calc_ringing(img, qmi );
   ring_part =  img .* ( (img<0) & ~qmi);
   rr = -sum( ring_part(:) )/sum( img(:).*qmi(:) );

function [xmean,ymean,equiv_circ,map,qmi,img] = calc_cofg(img);
%  if abs(max(img(:))) < abs(min(img(:))); img= -img; end
   qmi = calc_hm_set( img, 0.25 );
   if sum(img(:) & qmi(:))<0 ; keyboard ; end
   [x,y]=meshgrid(linspace(-1,1,32),linspace(-1,1,32)); map = x.^2+y.^2<1.1;
   qmi = qmi.*map; img = img.*map;

   ss_qmi = sum(qmi(:));
   xmean =  sum(sum( (qmi.*x) ))/ss_qmi; % centre of gravity
   ymean =  sum(sum( (qmi.*y) ))/ss_qmi;
   equiv_circ = (x-xmean).^2 + (y-ymean).^2 < ss_qmi/pi/(32/2)^2;
