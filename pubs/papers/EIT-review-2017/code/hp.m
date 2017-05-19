% Compare 2D algorithms
% $Id$
rng(0); % set seed for repeatable results

imb=  mk_common_model('c2C',16);

e= size(imb.fwd_model.elems,1);
bkgnd= 1;

% Solve Homogeneous model
img= mk_image(imb.fwd_model, bkgnd);
vh= fwd_solve( img );

% Add Two triangular elements
img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])=bkgnd * 2;
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=bkgnd * 2;
vi= fwd_solve( img );

% Add -12dB SNR
vi_n= vi; 
nampl= std(vi.meas - vh.meas)*10^(-18/20);
vi_n.meas = vi.meas + nampl *randn(size(vi.meas));

img.calc_colours.cb_shrink_move = [0.8 0.8 0];
show_fem(img,1);
axis square; axis off
print_convert('hp1.png')

% Compare 2D algorithms
% $Id$

for jj = [2 3 4];
  clear imgn
  switch(jj)
  case 4 % tik    noser lapl  nf   tv
    hpv = [0.05   0.25  0.03   0.5 1e-4];
  case 3
    hpv = [0.019  0.10  0.014  1.0 1e-5];
  case 2
    hpv = [0.0040 0.024 0.0028 3.0 3e-6];
  otherwise
    error('oops');
  end

  % Create Inverse Model
  inv2d= eidors_obj('inv_model', 'EIT inverse');
  inv2d.reconst_type= 'difference';
  inv2d.jacobian_bkgnd.value= 1;
  
  % This is not an inverse crime; inv_mdl != fwd_mdl
  imb=  mk_common_model('b2C',16);
  inv2d.fwd_model= imb.fwd_model;
  
  % Gauss-Newton solvers
  inv2d.solve=       @inv_solve_diff_GN_one_step;
  
  % Tikhonov prior
  inv2d.hyperparameter.value = hpv(1);
  inv2d.RtR_prior=   @prior_tikhonov;
  imgn(1)= inv_solve( inv2d, vh, vi_n);
  inva{1} = inv2d;
  
  % NOSER prior
  inv2d.hyperparameter.value = hpv(2);
  inv2d.RtR_prior=   @prior_noser;
  imgn(2)= inv_solve( inv2d, vh, vi_n);
  inva{2} = inv2d;
  
  % Laplace image prior
  inv2d.hyperparameter.value = hpv(3);
  inv2d.RtR_prior=   @prior_laplace;
  imgn(3)= inv_solve( inv2d, vh, vi_n);
  inva{3} = inv2d;
  
  % Automatic hyperparameter selection
  inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'value');
  inv2d.hyperparameter.func = @choose_noise_figure;
  inv2d.hyperparameter.noise_figure= hpv(4);
  inv2d.hyperparameter.tgt_elems= 1:4;
  inv2d.RtR_prior=   @prior_gaussian_HPF;
  inv2d.solve=       @inv_solve_diff_GN_one_step;
  imgn(4)= inv_solve( inv2d, vh, vi_n);
  inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'func');
  hp4 = choose_noise_figure(inv2d);
  inv2d.hyperparameter.value = hp4;
  inva{4} = inv2d;
  
  
  % Total variation using PDIPM
  inv2d.hyperparameter.value = hpv(5);
  inv2d.solve=       @inv_solve_TV_pdipm;
  inv2d.R_prior=     @prior_TV;
  inv2d.parameters.max_iterations= 10;
  inv2d.parameters.term_tolerance= 1e-3;
  
  %Vector of structs, all structs must have exact same (a) fields (b) ordering
  imgn5= inv_solve( inv2d, vh, vi_n);
  imgn5=rmfield(imgn5,'type'); imgn5.type='image';
  imgn(5)=imgn5;
  inva{5} = inv2d;
  
  % Output image
  for ii = 1:length(imgn)
    imgn(ii).calc_colours.ref_level = 0;
    imgn(ii).calc_colours.clim = 0.5;
    clf;
    h=show_fem(imgn(ii));
    set(h,'EdgeColor','none'); axis off;
    hold on; rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]); hold off;
%   print('-dpng',sprintf('hp%d%d.png',jj,ii));
    print_convert(sprintf('hp%d%d.png',jj,ii));

  end
  for ii = 1:length(inva)
     inv2d = inva{ii};
     inv2d.hyperparameter.tgt_elems = 1:4;
     [NF(ii,jj),SE(ii,jj)] =  calc_noise_figure(inv2d);
  end
     
end
clf; c=eidors_colourbar(imgn(1)); axis off;
set(c,'Location','SouthOutside','Position',[0.1143 0.1730 0.6750 0.0335]);
%print('-dpng','hp-cb.png');
print_convert('hp-cb.png');
