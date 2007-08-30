% Basic Image reconstruction
% $Id: tutorial110a.m,v 1.3 2007-08-30 03:58:27 aadler Exp $

% Load some data
load iirc_data_2006

% grey background
calc_colours('greylev',-.1);

% Get a 2D image reconstruction model
imdl= mk_common_model('c2c');
vi= real(v_rotate(:,9))/1e4; vh= real(v_reference)/1e4;

for idx= 1:3
   if     idx==1
      imdl.hyperparameter.value= 1e-4;
   elseif idx==2
      imdl.hyperparameter.value= 3e-3;
   elseif idx==3
      imdl.hyperparameter.value= 1e-1;
   end
   img= inv_solve(imdl, vh, vi);

   subplot(2,3,idx);
   show_slices(img);

   subplot(2,3,idx+3);
   z=calc_slices(img);
   c=calc_colours(z);
   mesh(z,c); view(-11,44);
 %NOT ONLY DO WE SPECIFY COLOURS, WE ACTUALLY WANT MATLAB TO USE THEM!
   set(gca,'Clim',[1 256]);
   set(gca,{'XLim','YLim','ZLim','XTickLabel','YTickLabel'}, ...
        {[1 64],[1 64],[-1.3,0.2],[],[]})
end


print -r100 -dpng tutorial110a.png;
