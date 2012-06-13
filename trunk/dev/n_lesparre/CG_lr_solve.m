function img= CG_lr_solve( inv_model, data)
% CG_LR_SOLVE absolute solver using the conjugate gradient method with
% L-Curve criterion. Convert the element conductivity into the
% logarithm of the corresponding resistivity

% img= CG_lr_solve( inv_model, data)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data      => EIT data
%
% The conjugate gradient method is iterative and computes at each iteration
% k the parameters:
% p(k+1) = p(k) + alpha(k+1)* delta(k+1)
% where alpha is the step length and delta the perturbation direction

% Parameters:
%     inv_model.parameters.max_iterations = N_max iter
%     inv_model.parameters.perturb = vector with at least 3 variables to
%     determine the step length
%     inv_model.parameters.lambda = vector in logarithm space to estimate 
%     regularization defined by the L-curve criterion
%     inv_model.parameters.normalisation = normalisation of the data sets,
%     for instance the geometrical factor used in the apparent resistivity
%     estimation

% (C) 2012 Nolwenn Lesparre. License: GPL version 2 or version 3
% $Id: CG_lr_solve.m 2596 2012-06-06 14:18:06Z nlesparre $

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

img = calc_jacobian_bkgnd( inv_model );

% Possibility to use a different parameterisation between the inverse and
% the forward problems. In that case, a coarse2fine matrix rely the
% elements to reconstruct the medium conductivity
if isfield(inv_model.fwd_model,'coarse2fine')
    nc = size(img.fwd_model.coarse2fine,2);
    img.elem_data = mean(img.elem_data)*ones(nc,1);
end

if isfield(inv_model.parameters,'max_iterations')
    iters = inv_model.parameters.max_iterations;
else
    iters = 1;
end

if ~isfield(inv_model.parameters,'perturb')
    img.parameters.perturb= [0.0001 0.001 0.01];
end

if ~isfield(inv_model.parameters,'lambda')
    img.parameters.lambda= logspace(-5,0,500);
end

if ~isfield(inv_model.parameters,'normalisation')
    img.parameters.normalisation= 1;
end

for k= 1:iters
    % Estimate residuals between data and estimation
    vsim=  fwd_solve(img);
    residuals = vsim.meas-data;
    residuals= inv_model.parameters.normalisation*residuals;
    
    % Calculate Jacobian
    disp(['Begin Jacobian computation - Iteration ' num2str(k)]);
    J = calc_jacobian( img );
    
    % Convert Jacobian as the adjusted parameters are the logarithm of the
    % resistivity
    img.logRes= -log(img.elem_data);
    dCond_dlogRes= -exp(-img.logRes);
    J = J.*repmat((dCond_dlogRes),1,size(data,1))';
    % Normalize the Jacobian
    J= inv_model.parameters.normalisation*J;
    
    % Some of the parameters may not be adjusted, as for instance the
    % background for huge forward models
    if isfield(inv_model.parameters,'fixed_background') && inv_model.parameters.fixed_background==1
        if isfield(inv_model.fwd_model.misc,'wall')
            J= J(:,1:nc-2);
        else
            J= J(:,1:nc-1);
        end
    end
    
    % Estimate the perturbation direction
    tol= svdAnalysisLcurvecrit(data,J);
    delta_params= pinv(J,tol)*residuals;
    
    if isfield(inv_model.parameters,'fixed_background') && inv_model.parameters.fixed_background==1
        if isfield(inv_model.fwd_model.misc,'wall')
            delta_params(nc-1:nc)= 0;
        else
            delta_params(nc)= 0;
        end
    end
    % Compute the step length and adjust the parameters
    img = line_optimize(img, delta_params, data);
end

end


function  img = line_optimize(imgk, dx, data1)
img = imgk;
perturb= img.fwd_model.parameters.perturb;

% Compute the forward model for eah perturbation step
mlist= perturb*0;
for i = 1:length(perturb);
    img.params_mapping.params= imgk.params_mapping.params + perturb(i)*dx;
    vsim = fwd_solve( img );
    dv = vsim.meas-data1;
    dv= img.parameters.normalisation*dv;
    mlist(i) = norm(dv);
end
% Select the best fitting step
pf = polyfit(perturb, mlist, 2);
fmin = -pf(2)/pf(1)/2; % poly minimum
fmin(fmin<0)= min(perturb);
if pf(1)*fmin^2+pf(2)*fmin+pf(3)>min(mlist) || fmin>max(perturb)
    [mk,ik]= min(cp);
    fmin= p(ik);
end
img.fwd_model.perturb= [fmin/4 fmin/2 fmin fmin*2 fmin*4] ;

% figure; semilogx(perturb,mlist,'xk',fmin,pf(1)*fmin^2+pf(2)*fmin+pf(3),'or'); hold on
% semilogx(p,pf(1)*p.^2+pf(2)*p+pf(3),'k','linewidth',2); axis tight
% xlabel('alpha','fontsize',30,'fontname','Times')
% ylabel('Normalized residuals','fontsize',30,'fontname','Times')
% title(['Best alpha = ' num2str(fmin,'%1.2e')],'fontsize',30,'fontname','Times')
% set(gca,'fontsize',30,'fontname','Times'); drawnow;

% Record the corresponding parameters
img.logRes= imgk.logRes + fmin*dx;
img.elem_data= exp(-img.logRes);
img.res_data= exp(img.logRes);
end


function tol= svdAnalysisLcurvecrit(data,img,J)
% Compute the singular value decompisition  
[U,S,V]= svd(J); svj= diag(S);
svj = svj(svj>eps);

% Estimate the model parameters from a forward model linear approximation
beta= U'*data;     %  adjust data
beta= beta(1:length(svj));

% Estimate the parameter of regularization
lambda= img.fwd_model.parameters.lambda;
% lambda= logspace(-4,-1,500);
% lambda= logspace(-6,-2,500);

% Determine the pinv tolerance following the L-curve criterion
SVJ= repmat(svj,1,length(lambda));
LAMBDA= repmat(lambda,length(svj),1);
FI= SVJ.^2./(SVJ.^2+LAMBDA.^2);
BETA= repmat(beta,1,length(lambda));
XL2= sum((FI.*BETA./SVJ).^2,1);
RES2= sum(((1-FI).*BETA).^2,1);
n= XL2; p= RES2;
nnp= (1-FI).*(FI.^2).*(BETA.^2)./SVJ.^2;
np= -(4./lambda).*sum(nnp,1);
kapa= (2*n.*p./np).*(lambda.^2.*np.*p+2*lambda.*n.*p+lambda.^4.*n.*np)./(lambda.^2.*n.^2+p.^2).^(3/2);

resi= sqrt(RES2); xlambda = sqrt(XL2);
[mk,ik]= min(kapa);
[m,ist]= min(abs(svj-lambda(ik)));
tol=svj(ist);

% figure;
% loglog(resi,xlambda,'k',resi(ik),xlambda(ik),'or','linewidth',2); axis tight; hold on
% yl= get(gca,'ylim');
% ylabel('Roughness','fontsize',30,'fontname','Times')
% xlabel('Residuals','fontsize',30,'fontname','Times')
% ylim(yl);
% title(['Best solution lambda=' num2str(lambda(ik),'%1.2e')],'fontsize',30,'fontname','Times')
% set(gca,'fontsize',30,'fontname','Times'); drawnow;
end  

function do_unit_test
unzip('../../htdocs/data_contrib/dg_geophysical_EIT/Mine_20FEV2004.zip')
data= load('Mine_20FEV2004_LI.tomel');
gps = load('Mine_20FEV2004.gps');
% Forward Model
 shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
              'solid mainobj= top and orthobrick(-100,-200,-100;425,10,100) -maxh=20.0;\n'];
 elec_pos = gps(:,2:4); e0 = elec_pos(:,1)*0;
 elec_pos = [  elec_pos, e0, e0+1, e0 ]; 
 elec_shape=[0.5,.5,.5];
 elec_obj = 'top';
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

% Load data and positions (unused in this tutorial)
  fmdl.stimulation = stim_meas_list( data(:,3:6) - 40100);

show_fem(fmdl);
%Reconstruction model
[cmdl]= mk_grid_model([], 2.5+[-50,-20,0:10:320,340,370], ...
                             -[0:2.5:10, 15:5:25,30:10:80,100,120]);
c2f = mk_coarse_fine_mapping( fmdl, cmdl);

fmdl.coarse2fine = c2f;
imdl= eidors_obj('inv_model','test');
imdl.fwd_model= fmdl;
imdl.rec_model= cmdl;
   imdl.reconst_type = 'difference';
   imdl.RtR_prior = @prior_laplace;
   imdl.solve = @inv_solve_diff_GN_one_step;
   imdl.hyperparameter.value = 0.1;
   imdl.fwd_model.normalize_measurements = 1;
   imdl.jacobian_bkgnd.value = 0.03;
% Difference image vs simulated data
vr = data(:,9);
img = mk_image( imdl );
vh = fwd_solve(img); vh = vh.meas;

imgr = inv_solve(imdl, vh, vr);

show_fem(imgr);
axis equal; axis tight
end
