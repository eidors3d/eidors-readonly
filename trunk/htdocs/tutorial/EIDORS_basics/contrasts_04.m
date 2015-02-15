sz = 5; img_idx = 'b';
for ellipse_x = [0.5,1,2];
   img = contrasts_04_modeller( sz, ellipse_x); 
   targ = img.fwd_model.mat_idx{1};
   for contrast = linspace( -3,3,7);
      img.elem_data( targ ) = exp(contrast);
      vv = fwd_solve(img);
      img_name = sprintf('04%c',img_idx); img_idx= img_idx+1;
      contrasts_03;
   end
end

