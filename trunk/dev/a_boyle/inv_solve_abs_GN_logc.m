function img= inv_solve_abs_GN_logc( inv_model, data1);
% INV_SOLVE_ABS_GN_LOGC
% This function calls INV_SOLVE_ABS_CORE to find a Gauss-Newton
% iterative solution using log conductivity.
%
% img = inv_solve_abs_GN_logc( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => EIT measurements
%
% See INV_SOLVE_ABS_CORE for arguments, options and parameters.
%
% (C) 2013 Alistair Boyle, Nolwenn Lespare, Andy Adler
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

% we assume resistivity output...
% probably should be conductivity to be consistent with the rest of EIDORS
%if ~isfield(inv_model.parameters, 'elem_output')
%   inv_model.parameters.elem_output = 'resistivity';
%end

% fixed working data... otherwise we wouldn't be calling this function!
inv_model.parameters.elem_working = 'log_conductivity';
img = inv_solve_abs_core(inv_model, data1);

function do_unit_test
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;410,10,100) -maxh=20.0;\n'];
e0 = linspace(0,310,64)';
elec_pos = [e0,0*e0,0*e0,1+0*e0,0*e0,0*e0];
elec_shape= [0.1,0.1,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
%fmdl.nodes = fmdl.nodes(:,[1,3,2]);

% spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17];
% multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1];
% fmdl.stimulation= stim_pattern_geophys( 64, 'Schlumberger', {'spacings', spacing,'multiples',multiples});

fmdl.stimulation= stim_pattern_geophys( 64, 'Wenner', {'spacings', 1:32} );

cmdl= mk_grid_model([], 2.5+[-30,5,20,30:10:290,300,315,340], ...
                            -[0:5:10 17 30 50 75 100]);
% cmdl = mk_grid_model([], 2.5+[-50,-20,0:10:310,330,360], ...
%                              -[0:2.5:10, 15:5:25,30:10:80,100,120]);
c2f = mk_coarse_fine_mapping( fmdl, cmdl);
S= sum(c2f,2); b= find(S<0.9999); a= find(S>=0.9999 & S<1);
c2f(a,:)= c2f(a,:)./repmat(S(a),1,size(c2f,2));
c2f(b,:)= 0; c2f(b,end+1)= 1;
fmdl.coarse2fine= c2f;


img = mk_image(fmdl,1);
fm_pts = interp_mesh(fmdl);
x_bary= fm_pts(:,1); z_bary= fm_pts(:,2);

z_params= (min(fmdl.nodes(:,2)):max(fmdl.nodes(:,2)))';
a = 0.36;
b = 130;
x_params= a*z_params+b;
xlim=interp1(z_params,x_params,z_bary);
img.elem_data(x_bary>xlim)= 0.01;

% img2= mk_image(fmdl,img.elem_data);
% figure; show_fem(img2);

% img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[100;-30;0;50])*100);
% img.elem_data(img.elem_data==0)= 0.1;
dd  = fwd_solve(img);
figure; show_fem(img); title('model');

imdl= eidors_obj('inv_model','test');
imdl.fwd_model= fmdl;
imdl.rec_model= cmdl;
imdl.fwd_model.normalize_measurements = 0;
imdl.rec_model.normalize_measurements = 0;
imdl.RtR_prior = @prior_laplace;
imdl.solve = @inv_solve_abs_GN_logc;
imdl.reconst_type = 'absolute';
imdl.hyperparameter.value = 0.1;
imdl.jacobian_bkgnd.value = 1;


img1= mk_image(fmdl,1);
vh1= fwd_solve(img1);
normalisation= 1./vh1.meas;
I= speye(length(normalisation));
I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

imdl.inv_solve.calc_solution_error = 0;

imdl.parameters.verbose = 3;
imdl.parameters.normalisation= I;
imdl.parameters.homogeneization= 1;
imdl.parameters.fixed_background= 1;
%imdl.parameters.perturb= [0 5*logspace(-5,-1,5)];
%imdl.parameters.max_iterations= 10;
%imdl.parameters.plot_line_optimize = 1;

imgr= inv_solve(imdl, dd);

imgGNd= imgr;
imgGNd.fwd_model.coarse2fine= cmdl.coarse2fine;
imgGNd.elem_data= log10(imgGNd.res_data(1:end-1));
imgGNd.calc_colours.clim= 1.5;
imgGNd.calc_colours.ref_level= 1.5;

elec_posn= zeros(length(fmdl.electrode),3);
for i=1:length(fmdl.electrode)
    elec_posn(i,:)= mean(fmdl.nodes(fmdl.electrode(1,i).nodes,:),1);
end

figure; show_fem(imgGNd,1);
hold on; plot(elec_posn(:,1),elec_posn(:,3),'k*');
axis tight; ylim([-100 0.5])
xlabel('X (m)','fontsize',20,'fontname','Times')
ylabel('Z (m)','fontsize',20,'fontname','Times')
set(gca,'fontsize',20,'fontname','Times');

img = mk_image( imdl );
img.elem_data= imgr.elem_data;
vCG= fwd_solve(img); vCG = vCG.meas;

figure; plot(I*(dd.meas-vCG)); title('data misfit');
figure; hist(abs(I*(dd.meas-vCG)),50); title('|data misfit|, histogram'); xlabel('|misfit|'); ylabel('count');

figure; show_pseudosection( fmdl, I*dd.meas); title('measurement data');
figure; show_pseudosection( fmdl, I*vCG); title('reconstruction data');
figure; show_pseudosection( fmdl, (vCG-dd.meas)./dd.meas*100); title('data misfit');
