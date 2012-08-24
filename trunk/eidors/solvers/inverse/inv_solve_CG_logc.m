function img= inv_solve_CG_logc( inv_model, data)
% inv_solve_CG_logc absolute solver using the conjugate gradient method with
% L-Curve criterion. Convert the element conductivity into their
% logarithm

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
% $Id$


% Possibility to use a different parameterisation between the inverse and
% the forward problems. In that case, a coarse2fine matrix rely the
% elements to reconstruct the medium conductivity

img = calc_jacobian_bkgnd( inv_model );
img.parameters= inv_model.parameters;

if isfield(inv_model.fwd_model,'coarse2fine')
    nc = size(img.fwd_model.coarse2fine,2);
    img.elem_data = mean(img.elem_data)*ones(nc,1);
end

if isfield(img.parameters,'max_iterations')
    iters = inv_model.parameters.max_iterations;
else
    iters = 1;
end

if ~isfield(img.parameters,'perturb')
    img.parameters.perturb= [0.0001 0.001 0.01];
end

if ~isfield(img.parameters,'lambda')
    img.parameters.lambda= logspace(-5,0,500);
end

if ~isfield(img.parameters,'normalisation')
    img.parameters.normalisation= 1;
end

if isfield(img.parameters,'homogeneization')
    img = homogeneous_estimate( img, data );
end



residuals= zeros(size(data,1),iters+1);
for k= 1:iters
    % Calculate Jacobian
    disp(['Begin Jacobian computation - Iteration ' num2str(k)]);
    J = calc_jacobian( img ); 

    % Convert Jacobian as the adjusted parameters are the logarithm of the
    % conductivity
    img.logCond= log(img.elem_data);
    dCond_dlogCond= img.elem_data;
    J = J.*repmat((dCond_dlogCond),1,size(data,1))';
    
    % Normalize the Jacobian
    J= img.parameters.normalisation*J;
       
    % Some of the parameters may not be adjusted, as for instance the
    % background for huge forward models
    if isfield(img.parameters,'fixed_background') && img.parameters.fixed_background==1
        if isfield(img.fwd_model.misc,'wall')
            J= J(:,1:nc-2);
        else
            J= J(:,1:nc-1);
        end
    end

    if k==1
        if isfield(img.parameters,'tol')
            tol= img.parameters.tol;
        else
            tol= svdAnalysisLcurvecrit(img.parameters.normalisation*data,img,J);
        end
    end
    % Estimate residuals between data and estimation
    vsim=  fwd_solve(img);
    res = data-vsim.meas;
    residuals(:,k)= img.parameters.normalisation*res;
    
    % Compute the step length for the Conjugate Gradient increment
    delta_params= pinv(J,tol)*(img.parameters.normalisation*(res));

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

vsim=  fwd_solve(img);
residuals(:,k+1) = img.parameters.normalisation*(vsim.meas-data);
img.residuals= residuals;
img.estimation= vsim.meas;

end


function  img = line_optimize(imgk, dx, data1)
img = imgk;
perturb= img.parameters.perturb;

% Compute the forward model for eah perturbation step
mlist= perturb*0;
for i = 1:length(perturb); 
    img.logCond= imgk.logCond + perturb(i)*dx;
    img.elem_data= exp(img.logCond);
    vsim = fwd_solve( img );
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

img.logCond= imgk.logCond + fmin*dx;
img.elem_data= exp(img.logCond);
vsim = fwd_solve( img );
dv = vsim.meas-data1;
dv= img.parameters.normalisation*dv;
% norm(dv)

if norm(dv) > mlist(1)
    img.parameters.perturb= [0 logspace(-4,-2,5)];
    img.logCond= imgk.logCond;
else
      img.parameters.perturb= [0 fmin/2 fmin fmin*2];
      img.logCond= imgk.logCond + fmin*dx;
      if size(imgk.fwd_model.elems,1) <= 200000
          img.parameters.perturb= [0 fmin/4 fmin/2 fmin fmin*2 fmin*4];
      else
          img.parameters.perturb= [0 fmin/2 fmin fmin*2];
      end
end

% figure; semilogx(perturb(2:end),mlist(2:end),'xk',fmin,pf(1)*log10(fmin)^2+pf(2)*log10(fmin)+pf(3),'or'); hold on
% semilogx(10.^p,pf(1)*(p).^2+pf(2)*p+pf(3),'k','linewidth',2); axis tight
% xlabel('alpha','fontsize',20,'fontname','Times')
% ylabel('Normalized residuals','fontsize',20,'fontname','Times')
% title({['Best alpha = ' num2str(fmin,'%1.2e')] ; ...
%     ['norm no move = ' num2str(mlist(1),4)]},'fontsize',30,'fontname','Times')
% set(gca,'fontsize',20,'fontname','Times'); drawnow;

% Record the corresponding parameters
img.elem_data= exp(img.logCond);
img.res_data= exp(-img.logCond);
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
% lambda= logspace(-4,-1,500);
% lambda= logspace(-6,-2,500);

% Determine the pinv tolerance following the L-curve criterion
SVJ= repmat(svj,1,length(lambda));
LAMBDA= repmat(lambda,length(svj),1);
FI= SVJ.^2./(SVJ.^2+LAMBDA.^2);
BETA= repmat(beta,1,length(lambda));
XL2= sum((FI.*BETA./SVJ).^2,1);
RES2= sum(((1-FI).*BETA).^2,1);
resi= sqrt(RES2); xlambda = sqrt(XL2);
% 
n= XL2; p= RES2;
nnp= (1-FI).*(FI.^2).*(BETA.^2)./SVJ.^2;
np= -(4./lambda).*sum(nnp,1);
kapa= (2*n.*p./np).*(lambda.^2.*np.*p+2*lambda.*n.*p+lambda.^4.*n.*np)./(lambda.^2.*n.^2+p.^2).^(3/2);

[mk,ik]= min(kapa);
[m,ist]= min(abs(svj-lambda(ik)));
tol=svj(ist);

n= length(lambda);

pf1 = polyfit(log10(resi(1:round(n*0.2))), log10(xlambda(1:round(n*0.2))), 1);
cp1= pf1(1)*log10(resi)+pf1(2);

pf2 = polyfit(log10(resi(round(n*0.8):end)), log10(xlambda(round(n*0.8):end)), 1);
cp2= pf2(1)*log10(resi)+pf2(2);
% 
[mk,ik]= min(abs(cp1-cp2));
[m,ist]= min(abs(svj-lambda(ik)));
tol= svj(ist);

% [mu,iu]= min(abs(img.parameters.tol-lambda));
% [ms,ist1]= min(abs(svj-1));
% % [ms,ist1]= min(abs(svj-1.5))
% [ms,ist2]= min(abs(svj-2));
% % [ms,ist1]= min(abs(svj-5))
% [ist1 ist2];

% figure;
% loglog(resi,xlambda,'k','linewidth',2); axis tight; hold on
% % [x,y]= ginput;
% % [mk,ik]= min(abs(resi-x));
% % [m,ist]= min(abs(svj-lambda(ik)));
% % tol= svj(ist)
% % tol= lambda(ik)
% % loglog(resi(ik),xlambda(ik),'or','linewidth',2)
% loglog(resi(iu),xlambda(iu),'xr','linewidth',2)
% 
% % loglog(resi,10.^cp1,resi,10.^cp2,'linewidth',2)
% yl= get(gca,'ylim');
% ylabel('Roughness','fontsize',20,'fontname','Times')
% xlabel('Residuals','fontsize',20,'fontname','Times')
% ylim(yl);
% % title({['Best solution lambda=' num2str(lambda(ik),'%1.2e')]; ...
% %     ['Number of Singular vector involved = ' num2str(ist)]},'fontsize',20,'fontname','Times')
% set(gca,'fontsize',20,'fontname','Times'); drawnow;
% im= 'Figures/LCurveGallery';
% saveas(gcf, [im '.pdf'], 'pdf');

% figure;
% semilogx(resi,kapa,'k','linewidth',2); axis tight; hold on
% yl= get(gca,'ylim');
% ylabel('Kapa','fontsize',30,'fontname','Times')
% xlabel('Residuals','fontsize',30,'fontname','Times')
% ylim(yl);
% title(['Best solution lambda=' num2str(lambda(ik),'%1.2e')],'fontsize',30,'fontname','Times')
% set(gca,'fontsize',30,'fontname','Times'); drawnow;


end  


function img = homogeneous_estimate( img, data )
%    img = calc_jacobian_bkgnd( imdl );
%    vs = fwd_solve(img);
% 
%    if isstruct(data)
%       data = data.meas;
%    else
%      meas_select = [];
%      try
%         meas_select = img.fwd_model.meas_select;
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
%       img.elem_data = mean(img.elem_data)*pf(1);
%    end
%     img.jacobian_bkgnd.value = mean(img.elem_data)*pf(1);
%     disp(['Homogeneous resistivity = ' num2str(1/(mean(img.elem_data)*pf(1)),3) ' Ohm.m'])

    conductivity_background= (img.parameters.normalisation*data)\(ones(size(data)));
    img.jacobian_bkgnd.value = conductivity_background;
    img.elem_data= ones(size(img.elem_data))*conductivity_background;
    disp(['Homogeneous resistivity = ' num2str(1/(conductivity_background),3) ' Ohm.m'])
end