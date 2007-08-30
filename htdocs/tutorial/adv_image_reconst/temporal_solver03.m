% Image reconstruction of moving objects $Id: temporal_solver03.m,v 1.1 2007-08-30 03:55:34 aadler Exp $

image_select= length(xyr_pt)/2+1;; % this image is at 9 O'Clock
time_steps=  3; ts_expand= 5;
time_weight= .8;
ts_vec= -time_steps:time_steps;

% vi_sel is the inhomog data used by the algorithm
% sel is the image to show (ie. last for kalman, middle for temporal)
for alg=1:4
   if     alg==1; % GN Solver
      im_sel = image_select;
      vi_sel = vi_n(:,im_sel);

      sel  = 1;
      imdl= imdl_GN;

   elseif alg==2  % Weighted GN Solver
      im_sel= image_select+ ts_vec*ts_expand;
      weight= (time_weight.^abs(ts_vec));
      vi_sel= vi_n(:,im_sel) * weight(:) / sum(weight);

      sel  = 1;
      imdl= imdl_GN;

   elseif alg==3  % Temporal Solver
      im_sel= image_select+ ts_vec*ts_expand;
      vi_sel= vi_n(:,im_sel);

      sel  = 1 + time_steps; % choose the middle
      imdl= imdl_TS;

   elseif alg==4  % Kalman Solver
      im_sel= image_select+ (-12:0)*ts_expand; %let Kalman warm up
      vi_sel= vi_n(:,im_sel);

      sel  = length(im_sel); % choose the last
      imdl= imdl_KS;

   end

   imdl.fwd_model.normalize_measurements= 0;
   img= inv_solve( imdl, vh, vi_sel);
   % only show image for sel (ie. last for kalman, middle for temporal)
   img.elem_data= img.elem_data(:,sel);

   show_fem(img);
   axis equal
   axis([-1.1, 0.1, -1.1, 1.1]);
   set(gca,{'XTicklabel','YTicklabel'},{'',''});
%
% Put circles where the data points used for reconstruction are
%
   theta= linspace(0,2*pi,length(xyr_pt)/ts_expand);
   xr=   cos(theta);       yr=   sin(theta);
   xpos= xyr_pt(1,im_sel); ypos= xyr_pt(2,im_sel);
   rad= xyr_pt(3,im_sel);
   hold on;
   for i=1:length(xpos)
       hh= plot(rad(i)*xr+ xpos(i),rad(i)*yr+ ypos(i));
       set(hh,'LineWidth',3,'Color',[0,0,0]);
   end
   hold off;

   print('-r75','-dpng',sprintf('temporal_solver03%c.png',96+alg));
end % for i
