function img= inv_solve_abs_CG_params( inv_model, data)
% INV_SOLVE_ABS_CG_PARAMS absolute solver using the conjugate gradient method with
% Least Square criterion. This function operates with a mapping function
% linking the inverse parameters to the forward elements. The mapping
% function is defined in the inv_model.params_mapping structure.

% img= CG_params_solve( inv_model, data)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data      => EIT data
%
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

% inv_model.params_mapping structure may at least contain:
% params_mapping.function
% params_mapping.params
% params_mapping.perturb

% (C) 2012 Nolwenn Lesparre. License: GPL version 2 or version 3
% $Id$


% Necessity to define a different parameterisation of the inverse problem
% by respect to the forward problem one. In that case, a mapping function
% rely the elements to reconstruct the medium conductivity

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

if isfield(inv_model,'params_mapping') && isfield(inv_model.params_mapping,'function')
    mapping_function= inv_model.params_mapping.function;
    img= feval(mapping_function,inv_model);
else
    error('The inverse model must contain a field "params_mapping" where a mapping function links the forward and inverse parameters');
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

np= size(img.params_mapping.params,1);
residuals= zeros(size(data,1),iters+1);
for k= 1:iters
    % Estimate residuals between data and estimation
    img= feval(mapping_function,img);
    vsim=  fwd_solve(img);
    res = data-vsim.meas;
    residuals(:,k)= img.parameters.normalisation*res;
%     norm(residuals(:,k))
    % Calculate jacobian
    disp(['Begin Jacobian computation - Iteration ' num2str(k)]);
%     eidors_cache clear jacobian
    J = calc_jacobian( img );
    % Normalize the Jacobian
    J= img.parameters.normalisation*J;
%     figure; plot(J); 
    if k==1;
    if np<=7
        [U,S,V]= svd(J); svj= diag(S);
        svj = svj(svj>eps);%figure; semilogy(svj)
%         tol= svj(np);
%         tol= svj(round(np/2));
        tol= svj(end)
    elseif np<=100
        [U,S,V]= svd(J); svj= diag(S);
        svj = svj(svj>eps);
        tol= svj(round(np/2));
    else
        tol= svdAnalysisLcurvecrit(img.parameters.normalisation*data,img,J);
%         tol= svdAnalysisLcurvecrit(residuals(:,k),img,J);
    end
    end
%     delta_params= pinv(J,tol)*residuals(:,k);
    delta_params= pinv(J,tol)*(img.parameters.normalisation*res);
    if k== 1; figure; plot(delta_params,'x-'); drawnow; end
    % Compute the step length and adjust the parameters
    img = line_optimize(img, delta_params, data);
    [(img.params_mapping.params(1)) ; img.params_mapping.params(2) ; exp(img.params_mapping.params(end-1:end))]
end
%     [(img.params_mapping.params(1)) ; img.params_mapping.params(2) ; ...
%         exp(img.params_mapping.params(end-1:end))]
% end
img.residuals= residuals;
end


function  img = line_optimize(imgk, dx, data1)
img = imgk;
perturb= img.parameters.perturb;
mapping_function= img.params_mapping.function;
% Compute the forward model for eah perturbation step
mlist= perturb*0;
for i = 1:length(perturb);
    img.params_mapping.params= imgk.params_mapping.params + perturb(i)*dx;
    img= feval(mapping_function,img);
    vsim = fwd_solve( img );
    dv = calc_difference_data( vsim , data1, img.fwd_model);% vsim.meas-data1;
    dv= img.parameters.normalisation*dv;
    mlist(i) = norm(dv);
end
% Select the best fitting step
pf= polyfit(log10(perturb(2:end)), mlist(2:end), 2);
fmin= 10.^(-pf(2)/pf(1)/2); % poly minimum
p= linspace(log10(perturb(2)),log10(perturb(end)),50);
cp= pf(1)*p.^2+pf(2)*p+pf(3);

if pf(1)*(log10(fmin))^2+pf(2)*log10(fmin)+pf(3) > min(mlist(2:end)) || log10(fmin)>max(perturb)
    [mk,ik]= min(cp);
    fmin= 10.^(p(ik));
end

img.params_mapping.params= imgk.params_mapping.params + perturb(i)*dx;
img= feval(mapping_function,img);
vsim = fwd_solve( img );
dv = vsim.meas-data1;
dv= img.parameters.normalisation*dv;
norm(dv)

% if pf(1)*(log10(fmin))^2+pf(2)*log10(fmin)+pf(3) > mlist(1)
if norm(dv) > mlist(1)
    img.parameters.perturb= [0 logspace(-4,-2,5)];
    img.params_mapping.params= imgk.params_mapping.params;
else
%       img.parameters.perturb= [0 fmin/2 fmin fmin*2];
    img.params_mapping.params= imgk.params_mapping.params + perturb(i)*dx;
    img.parameters.perturb= [0 fmin/4 fmin/2 fmin fmin*2 fmin*4];
end

img= feval(mapping_function,img);


figure;
semilogx(perturb,mlist,'xk',fmin,pf(1)*log10(fmin)^2+pf(2)*log10(fmin)+pf(3),'or'); hold on
semilogx(10.^p,pf(1)*(p).^2+pf(2)*p+pf(3),'k','linewidth',2); axis tight
% semilogx(perturb,mlist,'xk',fmin,pf(1)*fmin^2+pf(2)*fmin+pf(3),'or'); hold on
% semilogx(p,pf(1)*p.^2+pf(2)*p+pf(3),'k','linewidth',2); axis tight
xlabel('alpha','fontsize',30,'fontname','Times')
ylabel('Normalized residuals','fontsize',30,'fontname','Times')
title({['Best alpha = ' num2str(fmin,'%1.2e')] ; ...
    ['norm no move = ' num2str(round(mlist(1)))]},'fontsize',30,'fontname','Times')
set(gca,'fontsize',30,'fontname','Times'); drawnow;

end


function tol= svdAnalysisLcurvecrit(data,img,J)
% Compute the singular value decompisition
[U,S,V]= svd(J); svj= diag(S);
svj = svj(svj>eps);

% Estimate the model parameters from a forward model linear approximation
beta= U'*data;     %  adjust data
beta= beta(1:length(svj));

% Estimate the parameter of regularization
lambda= img.parameters.lambda;

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

figure;
loglog(resi,xlambda,'k',resi(ik),xlambda(ik),'or','linewidth',2); axis tight; hold on
yl= get(gca,'ylim');
ylabel('Roughness','fontsize',30,'fontname','Times')
xlabel('Residuals','fontsize',30,'fontname','Times')
ylim(yl);
title(['Best solution lambda=' num2str(lambda(ik),'%1.2e')],'fontsize',30,'fontname','Times')
set(gca,'fontsize',30,'fontname','Times'); drawnow;
end   


function do_unit_test
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;410,10,100) -maxh=20.0;\n'];
e0 = linspace(0,310,64)';
elec_pos = [e0,0*e0,0*e0,1+0*e0,0*e0,0*e0];
elec_shape= [0.1,0.1,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
fmdl.jacobian= @jacobian_diff_map_func;
%fmdl.nodes = fmdl.nodes(:,[1,3,2]);

% spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17];
% multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1];
% fmdl.stimulation= stim_pattern_geophys( 64, 'Schlumberger', {'spacings', spacing,'multiples',multiples});

fmdl.stimulation= stim_pattern_geophys( 64, 'Wenner', {'spacings', 1:32} );
img= mk_image(fmdl,1/20);
fm_pts= interp_mesh(fmdl);
x_bary= fm_pts(:,1); z_bary= fm_pts(:,2);

z_params= (min(fmdl.nodes(:,2)):max(fmdl.nodes(:,2)))';
a = 0.36; b = 130;
x_params= a*z_params+b;
xlim=interp1(z_params,x_params,z_bary);
% xlim= 130;
img.elem_data(x_bary>xlim)= 1/120;

dd  = fwd_solve(img);
% figure;  show_fem(img);

img1= mk_image(fmdl,1);
vh1= fwd_solve(img1);
normalisation= 1./vh1.meas;
I= speye(length(normalisation));
I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

a = 0.3; b = 150;
res_params= log([10 100]');

z_params= (min(fmdl.nodes(:,2)):max(fmdl.nodes(:,2)))';
x_params= a*z_params+b;

imdl = eidors_obj('inv_model','test');
imdl.fwd_model= fmdl;
imdl.params_mapping.params= [a; b ; res_params];
imdl.params_mapping.perturb= [0.005; 2; 0.05 ; 0.5];
% imdl.params_mapping.params= [b; res_params];
% imdl.params_mapping.perturb= [2; 0.05; 0.5];

imdl.params_mapping.function = @border_mapping;
imdl.params_mapping.data.x_bary = x_bary;
imdl.params_mapping.data.z_bary = z_bary;
imdl.params_mapping.data.res_params = res_params;
imdl.params_mapping.data.x_params = x_params;
imdl.params_mapping.data.z_params = z_params;
imdl.params_mapping.data.a = a;
imdl.params_mapping.data.b = b;
imdl.reconst_type= 'absolute';
imdl.solve = @inv_solve_abs_CG_params;
imdl.normalize_measurements= 1;
imdl.parameters.normalisation= I;
imdl.parameters.lambda= logspace(-5,5,1000);
imdl.parameters.max_iterations = 10;
imdl.parameters.perturb= [0 logspace(-2,0,5)];
imdl.jacobian_bkgnd.value = 1;
imgr= inv_solve(imdl, dd);

img = mk_image( imdl );
img.elem_data= imgr.elem_data;
vCG= fwd_solve(img); vCG = vCG.meas;
figure; hist(I*(dd.meas-vCG),50)
show_pseudosection( fmdl, I*dd.meas, 'HorizontalDownward')
show_pseudosection( fmdl, I*vCG, 'HorizontalDownward')
show_pseudosection( fmdl, (vCG-dd.meas)./dd.meas*100','HorizontalDownward')

end


function img = border_mapping(img)
%% Function to be called to perform the mapping in the forward problem
z= img.params_mapping.data.z_params;
res= img.params_mapping.params(end-1:end);

a= img.params_mapping.params(1);
b= img.params_mapping.params(2);
% b= img.params_mapping.params(1);
% 
% a= img.params_mapping.data.a;
% b= img.params_mapping.params(1);
x= z*a+b; 


% x2= img.params_mapping.data.x_params2;
% z2= img.params_mapping.data.z_params2;

xi= img.params_mapping.data.x_bary;
zi= img.params_mapping.data.z_bary;
xlim=interp1(z,x,zi);
% zlim=interp1(x2,z2,xi);

% figure; plot(xi,zlim,'*')
% hold on, plot(xlim,zi,'x')


vi= zeros(size(img.fwd_model.elems,1),1) + res(1);
% vi(xi>b)= res(2);

vi(xi>xlim)= res(2);
% vi(zi>zlim&xi>xlim+15)= res(1);

img.elem_data= exp(-vi);

% img.params_mapping.data.x_params = x;
% img.params_mapping.params= [x(2:end) ;  res];

img.params_mapping.data.x_params = x;
% img.params_mapping.params= [b ;  res];
img.params_mapping.params= [a ; b ;  res];
% img.params_mapping.params= [b ;  res];

end
