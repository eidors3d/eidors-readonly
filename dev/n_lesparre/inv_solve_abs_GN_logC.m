function img= inv_solve_abs_GN_logC( inv_model, data1);
% INV_SOLVE_ABS_GNR absolute solver using Gauss Newton approximation
% img= gn_abs_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

% Step 1: fit to background
img = calc_jacobian_bkgnd( inv_model );
img.log_conductivity.elem_data = img.elem_data;
img = rmfield(img, 'elem_data');
% img = homogeneous_estimate( inv_model, data1 );
if isfield(inv_model.fwd_model,'coarse2fine')
    nc = size(img.fwd_model.coarse2fine,2);
    img.log_conductivity.elem_data = mean(img.log_conductivity.elem_data)*ones(nc,1);
end

% Load the paramaters for the inversion
if isfield(inv_model,'parameters')
    img.parameters= inv_model.parameters;
else img.parameters.default= [];
end

if ~isfield(img.parameters,'normalisation')
    img.parameters.normalisation= 1;
end

if ~isfield(img.parameters,'perturb')
    img.parameters.perturb= [0 0.0001 0.001 0.01];
end


if isfield(img.parameters,'homogeneization')
    img = homogeneous_estimate( img, data1 );
end


hp  = calc_hyperparameter( inv_model );
if isfield(img.parameters,'fixed_background') && img.parameters.fixed_background==1
    inv_model2= inv_model;
    inv_model2.fwd_model.coarse2fine= inv_model.fwd_model.coarse2fine(:,1:end-1);
    RtR = calc_RtR_prior( inv_model2 );
else
    RtR = calc_RtR_prior( inv_model );
end


W   = calc_meas_icov( inv_model );
hp2RtR= hp*RtR;

iters = 1;
try 
   iters = inv_model.parameters.max_iterations; 
end

img0 = img;
% img0.logCond= log(img0.elem_data);

residuals= zeros(size(data1,1),iters+1);

for k = 1:iters  
    img = physics_data_mapper(img);
    vsim=  fwd_solve(img);

    res = img.parameters.normalisation*(data1-vsim.meas);
    residuals(:,k)=res;
   
  % Calculate Jacobian
    disp(['Begin Jacobian computation - Iteration ' num2str(k)]);
    J = calc_jacobian( img ); 
    img = physics_data_mapper(img,1);
    % Convert Jacobian as the adjusted parameters are the logarithm of the
    % conductivity
%     img.logCond= log(img.elem_data);
%     dCond_dlogCond= img.elem_data;
%     J = J.*repmat((dCond_dlogCond),1,size(data1,1))';
    
    % Normalize the Jacobian
    J= img.parameters.normalisation*J;
    if isfield(img.parameters,'fixed_background') && img.parameters.fixed_background==1
        J= J(:,1:nc-1);
    end


        
    if isfield(img.parameters,'fixed_background') && img.parameters.fixed_background==1
        RDx = hp2RtR*(img0.log_conductivity.elem_data(1:end-1) ... 
                - img.log_conductivity.elem_data(1:end-1));
    else
        RDx = hp2RtR*(img0.log_conductivity.elem_data - img.log_conductivity.elem_data);
    end
%     figure; plot(RDx);
    dx = (J'*W*J + hp2RtR)\(J'*res + RDx);
    if isfield(inv_model.parameters,'fixed_background') && inv_model.parameters.fixed_background==1
        dx(nc)= 0;
    end
    
    img= line_optimize(img, dx, data1);
%         save new
%     return
%     img0= img;
end

img = physics_data_mapper(img);
vsim=  fwd_solve(img);
img = physics_data_mapper(img,1);
residuals(:,k+1) = img.parameters.normalisation*(vsim.meas-data1);
img.residuals= residuals;
img.estimation= vsim.meas;

% Fit a parabola to the linefit and pick the best point
% This is better than doing an exhaustive search
% function  img = line_optimize(imgk, dx, data1);
%   flist = [ 0.1,  0.5, 1.0];
%   clim = mean(imgk.elem_data)/10; % prevent zero and negative conductivity
%   img = imgk;
%   for i = 1:length(flist);
%      img.elem_data = imgk.elem_data + flist(i)*dx;
%      img.elem_data(img.elem_data <= clim ) = clim;
%      vsim = fwd_solve( img );
%      dv = calc_difference_data( vsim , data1, img.fwd_model);
%      mlist(i) = norm(dv);
%   end
%   pf = polyfit(flist, mlist, 2);
%   fmin = -pf(2)/pf(1)/2; % poly minimum
%   fmin(fmin>1) = 1; fmin(fmin<0) = 0;
% 
%   img.elem_data = imgk.elem_data + flist(i)*dx;
%   img.elem_data(img.elem_data <= clim ) = clim;
function  img = line_optimize(imgk, dx, data1)
img = imgk;
perturb= img.parameters.perturb;
% Compute the forward model for eah perturbation step
mlist= perturb*0;
for i = 1:length(perturb); 
    img.log_conductivity.elem_data= imgk.log_conductivity.elem_data + perturb(i)*dx;
    img = physics_data_mapper(img);
    vsim = fwd_solve( img );
    img = physics_data_mapper(img,1);
    dv = vsim.meas-data1;
    dv= img.parameters.normalisation*dv;
    mlist(i) = norm(dv);
end
% Select the best fitting step
pf= polyfit(log10(perturb(2:end)), mlist(2:end), 2);
fmin= 10.^(-pf(2)/pf(1)/2); % poly minimum
p= linspace(log10(perturb(2)),log10(perturb(end)),50);
cp= pf(1)*p.^2+pf(2)*p+pf(3);

if pf(1)*(log10(fmin))^2+pf(2)*log10(fmin)+pf(3) > min(mlist(2:end)) || log10(fmin)>max(perturb) || log10(fmin)<min(perturb(2:end)) 
    [mk,ik]= min(mlist(2:end));
    fmin= perturb(ik+1);
end

img.log_conductivity.elem_data = imgk.log_conductivity.elem_data + fmin*dx;
% img.elem_data= exp(img.logCond);
img = physics_data_mapper(img);
vsim = fwd_solve( img );
img = physics_data_mapper(img,1);
dv = vsim.meas-data1;
dv= img.parameters.normalisation*dv;

if norm(dv) > mlist(1)
    img.parameters.perturb= [0 logspace(-4,-2,5)];
    img.log_conductivity.elem_data= imgk.log_conductivity.elem_data;
else
      img.parameters.perturb= [0 fmin/2 fmin fmin*2];
      img.log_conductivity.elem_data= imgk.log_conductivity.elem_data + fmin*dx;
      if size(imgk.fwd_model.elems,1) <= 200000
          img.parameters.perturb= [0 fmin/4 fmin/2 fmin fmin*2 fmin*4];
      else
          img.parameters.perturb= [0 fmin/2 fmin fmin*2];
      end
end
save imgk
% figure; semilogx(perturb(2:end),mlist(2:end),'xk',fmin,pf(1)*log10(fmin)^2+pf(2)*log10(fmin)+pf(3),'or'); hold on
% semilogx(10.^p,pf(1)*(p).^2+pf(2)*p+pf(3),'k','linewidth',2); axis tight
% xlabel('alpha','fontsize',20,'fontname','Times')
% ylabel('Normalized residuals','fontsize',20,'fontname','Times')
% title({['Best alpha = ' num2str(fmin,'%1.2e')] ; ...
%     ['norm no move = ' num2str(mlist(1),4)]},'fontsize',30,'fontname','Times')
% set(gca,'fontsize',20,'fontname','Times'); drawnow;

% Record the corresponding parameters
% img.elem_data= exp(img.logCond);
% img.res_data= exp(-img.logCond);

function img = homogeneous_estimate( img, data )
%    img = calc_jacobian_bkgnd( imdl );
%    vs = fwd_solve(img);
% 
%    if isstruct(data)
%       data = data.meas;
%    else
%      meas_select = [];
%      try
%         meas_select = imdl.fwd_model.meas_select;
%      end
%      if length(data) == length(meas_select)
%         data = data(meas_select);
%      end
%    end
% 
%    pf = polyfit(data,vs.meas,1);
% 
%    if isfield(img.fwd_model,'coarse2fine');
% % TODO: the whole coarse2fine needs work here.
% %   what happens if c2f doesn't cover the whole region
%       nc = size(img.fwd_model.coarse2fine,2);
%       img.elem_data = mean(img.elem_data)*ones(nc,1)*pf(1);
%    else
%       img.elem_data = img.elem_data*pf(1);
%    end

 conductivity_background= (img.parameters.normalisation*data)\(ones(size(data)));
 img.homogenous_estimate = conductivity_background;
 img.log_conductivity.elem_data= ones(size(img.log_conductivity.elem_data,1),1) * log(conductivity_background);
 disp(['Homogeneous resistivity = ' num2str(1/(conductivity_background),3) ' Ohm.m'])
 
function do_unit_test
   unit_test_simdata
 
function unit_test_simdata
   [fmdl, cmdl] = unit_test_models;
   dd = unit_test_do_fwd_sim(fmdl);


   imdl= eidors_obj('inv_model','test');
   imdl.fwd_model= fmdl;
   imdl.rec_model= cmdl;
   imdl.RtR_prior = @prior_laplace;
   imdl.solve = @inv_solve_abs_GN_logC;
   imdl.reconst_type = 'absolute';
   imdl.hyperparameter.value = 0.1;
   imdl.jacobian_bkgnd.log_conductivity.elem_data = log(1);


   img1= mk_image(fmdl,1); vh1= fwd_solve(img1); v_m = length(vh1.meas);
   imdl.parameters.normalisation= spdiags(1./vh1.meas,0,v_m,v_m);
   imdl.parameters.homogeneization= 1;
   imdl.parameters.fixed_background= 1;
   imdl.parameters.perturb= [0,logspace(-5,-3,5)];
   imdl.parameters.max_iterations= 10;

   imdl.fwd_model.jacobian = @jacobian_log_conductivity;
   imdl.fwd_model.solve    = @fwd_solve_log_conductivity;

   imgr= inv_solve(imdl, dd);

   imgGNd= imgr;
   imgGNd.fwd_model.coarse2fine= cmdl.coarse2fine;
   % last value indicates background?
   imgGNd.elem_data= log10(exp(-imgGNd.log_conductivity.elem_data(1:end-1)));
   imgGNd.calc_colours.clim= 1.5;
   imgGNd.calc_colours.ref_level= 1.5;

      
   show_fem(imgGNd,1);

   elec_posn= zeros(length(fmdl.electrode),3);
   for i=1:length(fmdl.electrode)
       elec_posn(i,:)= mean(fmdl.nodes(fmdl.electrode(1,i).nodes,:),1);
   end
   hold on; plot(elec_posn(:,1),elec_posn(:,3),'k*');
   axis tight; ylim([-100 0.5])
   xlabel('X (m)','fontsize',20,'fontname','Times')
   ylabel('Z (m)','fontsize',20,'fontname','Times')
   set(gca,'fontsize',20,'fontname','Times');

   img = mk_image( imdl );
   img.log_conductivity.elem_data= imgr.log_conductivity.elem_data;
   img = physics_data_mapper(img);
   vCG= fwd_solve(img); 
   img = physics_data_mapper(img,1);
   vCG = vCG.meas;

   I= imdl.parameters.normalisation;
   figure; plot(I*(dd.meas-vCG))
   figure; hist(I*(dd.meas-vCG),50)

   show_pseudosection( fmdl, I*dd.meas, 'profile')
   show_pseudosection( fmdl, I*vCG, 'profile')
   show_pseudosection( fmdl, (vCG-dd.meas)./dd.meas*100, 'profile')


function [fmdl, cmdl] = unit_test_models
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

   fmdl.normalize_measurements = 0;
   cmdl.normalize_measurements = 0;

function dd = unit_test_do_fwd_sim(fmdl)
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
   % show_fem(img);
   %%
