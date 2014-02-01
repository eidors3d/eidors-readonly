function params = eval_GREIT_fig_merit(imgs, xyzr_pt, SimInLayer)
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
% SimInLayer: number of simulated points in one layer
% to simulate multi layers. Normalize AR value for the first layer

% (C) 2009 Andy Adler. Licensed under GPL v2 or v3
% $Id: eval_GREIT_fig_merit.m 3934 2013-04-20 15:46:00Z aadler $

% Zhao: AR was normalized to the targets at the center of the electrode plane

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
if N_imgs > 10 % doesn't make sense to normalize otherwise
    ctr_pts = sum((xyzr_pt(1:mdl_dim(mdl(1)),:)-repmat(ctr',1,size(xyzr_pt,2))).^2) < (0.05*r)^2;
    if any(ctr_pts)
        if nargin < 3
            params(1,:) = params(1,:)/mean(params(1,ctr_pts));
        else
            params(1,:) = params(1,:)/mean(params(1,ctr_pts(1:SimInLayer)));
        end
    else
        eidors_msg('eval_GREIT_fig_merit: no centre points found to normalize',1);
    end
end





function ampl = calc_amplitude(img)
   ampl = sum(img(:));
%    ampl = sum(img(img(:)>0)); % use only positive values // Zhao

function pe   = calc_posn_error(qmi, xmean, ymean, xy)
   % This definition allows + and - PE, but can also give zero in unexpected places
   pe = sqrt(sum(xy.^2)) - sqrt( xmean^2 + ymean^2);
   % This definition gives the absolute PE, but can't be negative
%  pe = sqrt((xy(1,:) - xmean).^2 + ...
%            (xy(2,:) - ymean).^2);

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
   if sum(img(:) & qmi(:))<0 ; 
     error('problem in CofG calculation');
   end
   
   pix_sz = (max(x(:)) - min(x(:))) *( max(y(:)) - min(y(:))) /numel(img);

   %map = x.^2+y.^2<1.1;
   qmi = qmi.*map; img = img.*map;

%  qmi = qmi .* img;  %USE THE IMAGE AMPLITUDE FOR THE CALCULATION

   ss_qmi = sum(qmi(:));
   xmean =  sum(sum( (qmi.*x) ))/ss_qmi; % centre of gravity
   ymean =  sum(sum( (qmi.*y) ))/ss_qmi;
   equiv_circ = (x-xmean).^2 + (y-ymean).^2 < pix_sz*ss_qmi/pi;
