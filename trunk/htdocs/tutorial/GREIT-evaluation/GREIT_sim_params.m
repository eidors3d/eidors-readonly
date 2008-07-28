function params = GREIT_sim_params(imgs, xyzr_pt)
% params = GREIT_sim_params(imgs)
%  params(1,:) = Image Amplitude
%  params(2,:) = Position Error => + toward centre, - toward edge
%  params(3,:) = Resolution

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

[x,y]=meshgrid(linspace(-1,1,32),linspace(-1,1,32)); map = x.^2+y.^2<1.1;

N_imgs = size(imgs,3);
for i= 1:N_imgs
   img= imgs(:,:,i); qmi = calc_hm_set( img, 0.25 );
   params(1,i) = calc_amplitude( img );
   params(2,i) = calc_posn_error( qmi, x,y, xyzr_pt(1:2,i) );
   params(3,i) = calc_resolution( qmi, map );
end

% TODO: Fix this when we start to care about units
params(1,:) = params(1,:)/mean(params(1,1:10));

function ampl = calc_amplitude(img)
   ampl = sum(img(:));

function pe   = calc_posn_error(qmi, x, y, xy)
   ss_qmi = sum(qmi(:));
   xmean = -sum(sum( (qmi.*x) ))/ss_qmi; % centre of gravity
   ymean = -sum(sum( (qmi.*y) ))/ss_qmi;

   pe = sqrt(sum(xy.^2)) - sqrt( xmean^2 + ymean^2);

function res  = calc_resolution(qmi, map)
   res = sqrt( sum(qmi(:)) / sum(map(:)));
