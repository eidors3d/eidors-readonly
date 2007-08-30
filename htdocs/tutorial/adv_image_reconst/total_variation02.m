% Compare 2D algorithms
% $Id: total_variation02.m,v 1.1 2007-08-30 03:30:50 aadler Exp $

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;

imb=  mk_common_model('c2c',16); %576 Elem model
inv2d.fwd_model= imb.fwd_model;
inv2d.fwd_model.misc.perm_sym= '{y}';

% Guass-Newton solvers
inv2d.solve=       @np_inv_solve;

% NOSER prior
inv2d.hyperparameter.value = 3e-1;
inv2d.RtR_prior=   @noser_image_prior;

imgr= inv_solve( inv2d, v_homg, v_simu);

%Simulation image
subplot(221)
show_slices(sim_img)
subplot(223)
z=calc_slices(sim_img);
c=calc_colours(z); mesh(z,c);
view(173,34);
set(gca,{'XLim','YLim','ZLim','XTickLabel','YTickLabel'}, ...
        {[1 64],[1 64],[0.9,1.1],[],[]})

%Reconstructed image
subplot(222)
show_slices(imgr)
subplot(224)
z=calc_slices(imgr);
c=calc_colours(z); mesh(z,c);
view(173,34);
set(gca,{'XLim','YLim','ZLim','XTickLabel','YTickLabel'}, ...
        {[1 64],[1 64],[-0.1,0.1],[],[]})

print -r100 -dpng total_variation02a.png;

