% TV Solutions % $Id: total_variation03.m,v 1.3 2007-08-30 03:58:28 aadler Exp $

% Create TV Inverse Model
invtv= eidors_obj('inv_model', 'EIT inverse');
invtv.reconst_type= 'difference';
invtv.jacobian_bkgnd.value= 1;

invtv.hyperparameter.value = 1e-3;
invtv.solve=       @ab_tv_diff_solve;
invtv.R_prior=     @ab_calc_tv_prior;
invtv.parameters.term_tolerance= 1e-3;
invtv.parameters.keep_iterations= 0;

invtv.fwd_model= inv2d.fwd_model;
   


maxiters= [1,3,6,15];
for i= 1:length(maxiters)
   invtv.parameters.max_iterations= maxiters(i);
   imgtv= inv_solve( invtv, v_homg, v_simu);

   %Reconstructed image
   subplot(2,length(maxiters),i);
   show_slices(imgtv)
   subplot(2,length(maxiters),i+length(maxiters));
   z=calc_slices(imgtv);
   c=calc_colours(z); mesh(z,c);
   title(sprintf('TV iters=%d',maxiters(i)))
   view(173,34);
   set(gca,{'XLim','YLim','ZLim','XTickLabel','YTickLabel'}, ...
           {[1 64],[1 64],[0.98,1.1],[],[]})
end

print -r150 -dpng total_variation03a.png;

