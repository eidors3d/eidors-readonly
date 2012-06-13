function img= CG_params_solve( inv_model, data)
% CG_PARAMS_SOLVE absolute solver using the conjugate gradient method with
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
% $Id: CG_lr_solve.m 2596 2012-06-06 14:18:06Z nlesparre $


% Necessity to define a different parameterisation of the inverse problem
% by respect to the forward problem one. In that case, a mapping function
% rely the elements to reconstruct the medium conductivity

mapping_function= inv_model.params_mapping.function;
img= feval(mapping_function,inv_model);

np= size(img.params_mapping.params,1);

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
    vh= vsim.meas;
    residuals =  calc_difference_data( vsim , data, img.fwd_model);%vsim.meas-data;
    residuals= inv_model.parameters.normalisation*residuals;

    % Calculate jacobian
    disp(['Begin Jacobian computation - Iteration ' num2str(k)]);
    J = calc_jacobian( img );
    
    % Normalize the Jacobian
    J= inv_model.parameters.normalisation*J;
    
    if np<=100
        [U,S,V]= svd(J); svj= diag(S);
        svj = svj(svj>eps);
        tol= svj(round(np/2));
    else
        tol= svdAnalysisLcurvecrit(residuals,img,J);
    end
    delta_params= pinv(J,tol)*residuals;
    % Compute the step length and adjust the parameters
    img = line_optimize(img, delta_params, data);
end

end


function  img = line_optimize(imgk, dx, data1)
img = imgk;
perturb= img.parameters.perturb;

% Compute the forward model for eah perturbation step
mlist= perturb*0;
for i = 1:length(perturb);
    img.params_mapping.params= imgk.params_mapping.params + perturb(i)*dx;
    vsim = fwd_solve( img );
    dv = calc_difference_data( vsim , data1, img.fwd_model);% vsim.meas-data1;
    dv= img.parameters.normalisation*dv;
    mlist(i) = norm(dv);
end
% Select the best fitting step
pf = polyfit(perturb, mlist, 2);
fmin = -pf(2)/pf(1)/2; % poly minimum
fmin(fmin<0)= min(perturb);
p= perturb(1):(perturb(end)-perturb(1))/30:perturb(end);
cp= pf(1)*p.^2+pf(2)*p+pf(3);

if pf(1)*fmin^2+pf(2)*fmin+pf(3)>min(mlist) || fmin>max(perturb)
    [mk,ik]= min(cp);
    fmin= p(ik);
end
if fmin>min(mlist); 
    [mk,ik]= min(mlist);
    fmin= perturb(ik);
end

img.parameters.perturb= [0 fmin/2 fmin fmin*2] ;
if fmin==0;
img.parameters.perturb= [0 logspace(-7,-1,7)] ;
end
% figure; semilogx(perturb,mlist,'xk',fmin,pf(1)*fmin^2+pf(2)*fmin+pf(3),'or'); hold on
% semilogx(p,pf(1)*p.^2+pf(2)*p+pf(3),'k','linewidth',2); axis tight
% xlabel('alpha','fontsize',30,'fontname','Times')
% ylabel('Normalized residuals','fontsize',30,'fontname','Times')
% title(['Best alpha = ' num2str(fmin,'%1.2e')],'fontsize',30,'fontname','Times')
% set(gca,'fontsize',30,'fontname','Times'); drawnow;

% Record the corresponding parameters
img.params_mapping.params= imgk.params_mapping.params + fmin*dx;
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

% figure;
% loglog(resi,xlambda,'k',resi(ik),xlambda(ik),'or','linewidth',2); axis tight; hold on
% yl= get(gca,'ylim');
% ylabel('Roughness','fontsize',30,'fontname','Times')
% xlabel('Residuals','fontsize',30,'fontname','Times')
% ylim(yl);
% title(['Best solution lambda=' num2str(lambda(ik),'%1.2e')],'fontsize',30,'fontname','Times')
% set(gca,'fontsize',30,'fontname','Times'); drawnow;
end   

