% TV Solutions % $Id$

% Create TV Inverse Model
invtv= eidors_obj('inv_model', 'EIT inverse');
invtv.reconst_type= 'difference';
invtv.jacobian_bkgnd.value= 1;

invtv.hyperparameter.value = 1e-3;
invtv.solve=       @TV_diffusivity_solve;
invtv.R_prior=     @calc_TV_prior;
invtv.parameters.term_tolerance= 1e-3;
invtv.parameters.keep_iterations= 1;

invtv.fwd_model= inv2d.fwd_model;
   
invtv.parameters.max_iterations= 20;
imgtv= inv_solve( invtv, v_homg, v_simu);

clf;
show_slices(imgtv)
print -r100 -dpng total_variation04a.png;

imgs= calc_slices(imgtv);
idx = [1,3,6,10,20];

subplot(211);
plot(squeeze(imgs(:,32,[idx])));
axis([1 64 -0.04 0.12]); set(gca,'XTickLabel',[]);
 legend('1','3','6','10','20');

subplot(212);
plot(squeeze(imgs(32,:,[idx])));
axis([1 64 -0.04 0.12]); set(gca,'XTickLabel',[]);

print -r75 -dpng total_variation04b.png;
