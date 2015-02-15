sz = 5;
for ellipse_x = [0.5,1,2];
   img = contrasts_04_modeller( sz, ellipse_x); 
   targ = img.fwd_model.mat_idx{1};
   for contrast = linspace( -2,2,5);
      img.elem_data( targ ) = exp(contrast);
      vh = fwd_solve(img);
      contrasts_03;
   end
end

