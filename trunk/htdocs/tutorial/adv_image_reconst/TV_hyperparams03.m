
alpha1list= [1:0.4:5];
alpha2list= [3:0.4:8];
for alpha2= alpha2list
   for alpha1= alpha1list
      invtv.hyperparameter.value =   10^-alpha2;
      invtv.ab_tv_diff_solve.alpha1= 10^-alpha1;
      imgs= inv_solve(invtv,vh,vi);
      imgs.elem_data= imgs.elem_data(:,[1,2,5,maxit]);
      imgs.calc_colours.window_range=.4;
      imgs.calc_colours.clim=1.4;
      out_img= show_slices(imgs);
      imwrite(out_img, colormap, sprintf(name_base,alpha1,alpha2) );
    end
end
 
