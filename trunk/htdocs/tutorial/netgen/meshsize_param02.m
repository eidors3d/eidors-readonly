% Reconstruct images $Id$
imdl= mk_common_model('c2c2',16);
imdl.RtR_prior= @laplace_image_prior;
imdl.hyperparameter.value = 0.01;

th= linspace(0,2*pi,100); th(end)= [];
xcirc= 0.5+0.1*cos(th); ycirc= 0.1*sin(th);

subplot(121);
show_fem(img); view(90,60);

subplot(122);
show_fem(inv_solve( imdl, vh, vi(1)));
line(xcirc,ycirc,'Color',[0,0.5,0],'LineWidth',2);
axis image
ylabel(sprintf('No. Elems= %d', size(img.fwd_model.elems,1)));

print -dpng -r100 meshsize_param02a.png
