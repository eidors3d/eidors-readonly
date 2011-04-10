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
mdl = imgs.fwd_model;
imgs = calc_slices(imgs);
map = ~isnan(squeeze(imgs(:,:,1))); %assume all imgs are the same shape
imgs(isnan(imgs)) = 0;
sz = size(imgs,1);
[x,y,bb_min,bb_max]=prepare_grid(sz,mdl);

N_imgs = size(imgs,3);
for i= 1:N_imgs
   [xmean,ymean,equiv_circ,qmi,img] = calc_cofg(imgs(:,:,i),map,x,y);
   params(1,i) = calc_amplitude( img );
   params(2,i) = calc_posn_error( qmi, xmean, ymean, xyzr_pt(1:2,i) );
   params(3,i) = calc_resolution( qmi, map );
   params(4,i) = calc_shape_deform( qmi, equiv_circ );
   params(5,i) = calc_ringing( img, qmi );
end

% TODO: Fix this when we start to care about units
ctr = bb_min + 0.5*(bb_max-bb_min);
r = max(0.5*(bb_max-bb_min));
ctr_pts = sum((xyzr_pt(1:2,:)-repmat(ctr',1,size(xyzr_pt,2))).^2) < (0.05*r)^2;
%params(1,:) = params(1,:)/mean(params(1,1:10));  
params(1,:) = params(1,:)/mean(params(1,ctr_pts));



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
   
function [x,y,bb_min,bb_max]=prepare_grid(sz,mdl)
      % bounding box
   bnd = unique(mdl.boundary);
   bb_min = min(mdl.nodes(bnd,:));
   bb_max = max(mdl.nodes(bnd,:));
   
   [x,y]=meshgrid(linspace(bb_min(1),bb_max(1),sz),linspace(bb_min(2),bb_max(2),sz)); 
   
function [xmean,ymean,equiv_circ,qmi,img] = calc_cofg(img,map,x,y);
%  if abs(max(img(:))) < abs(min(img(:))); img= -img; end
   sz = size(img,1);
   qmi = calc_hm_set( img, 0.25 );
   if sum(img(:) & qmi(:))<0 ; keyboard ; end
   
   

   %map = x.^2+y.^2<1.1;
   qmi = qmi.*map; img = img.*map;

   ss_qmi = sum(qmi(:));
   xmean =  sum(sum( (qmi.*x) ))/ss_qmi; % centre of gravity
   ymean =  sum(sum( (qmi.*y) ))/ss_qmi;
   equiv_circ = (x-xmean).^2 + (y-ymean).^2 < ss_qmi/pi/(sz/2)^2;
