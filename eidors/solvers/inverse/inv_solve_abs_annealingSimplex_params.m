function [img]= inv_solve_abs_annealingSimplex_params(inv_model, data)
% INV_SOLVE_ABS_ANNEALINGSIMPLEX_PARAMS absolute solver using the simplex annealing method. 
% This function operates with a mapping function linking the inverse
% parameters linking the inverse parameters to the forward elements.
% The mapping function is defined in the inv_model.params_mapping structure.

% img= annealingSimplex( inv_model, data)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data      => EIT data
%
% Parameters:
%     inv_model.parameters.tempInit = N_max iter
%     inv_model.parameters.tempfinale = vector with at least 3 variables
%     inv_model.parameters.cooldown = 
%     inv_model.parameters.normalisation

% (C) 2012 Nolwenn Lesparre. License: GPL version 2 or version 3
% $Id$


% Necessity to define a different parameterisation of the inverse problem
% by respect to the forward problem one. In that case, a mapping function
% rely the elements to reconstruct the medium conductivity

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

if isfield(inv_model,'params_mapping') &&  isfield(inv_model.params_mapping,'function')
    mapping_function= inv_model.params_mapping.function;
    img= feval(mapping_function,inv_model);
else
    error('The inverse model must contain a field "params_mapping" where a mapping function links the forward and inverse parameters');
end

if isfield(inv_model.parameters,'temp')
    tempInit = inv_model.parameters.temp;
    tempfinale= inv_model.parameters.tempfinale;
    cooldown= inv_model.parameters.cooldown;
    nMetro= inv_model.parameters.nMetro;
else
    tempInit = 1000;
    tempfinale= 0.01;
    cooldown= 0.95;
    nMetro= 3;
end

if ~isfield(inv_model.parameters,'normalisation')
    img.parameters.normalisation= 1;
end
temp= tempInit; k= 1;
while temp>=tempfinale
    temp= temp*cooldown;
    k=k+1;
end
D= zeros(k,1); temperature= zeros(k,1);
niter= k
temp= tempInit; k= 1;

% Estimate the number of parameters to adjust
np= size(img.params_mapping.params,1);

% Generates np+1 models with a Gaussian law with a mean equals to
% inv_model.params_mapping.params and a standard deviation of
% inv_model.params_mapping.perturb
if isfield(inv_model.params_mapping,'inital_model')
    modeles= inv_model.params_mapping.inital_model;
else
    modeles= randn(np,np+1).*repmat(inv_model.params_mapping.perturb,1,np+1) + ...
        repmat(inv_model.params_mapping.params,1,np+1);
    modelesr= modeles(end,:);
    keep= ismember(modeles(end,:)<modeles(end-1,:),1);
    modeles(end,keep)= modeles(end-1,keep);
    modeles(end-1,keep)= modelesr(keep);
end
modelesT= [];
costhhi= [];
modelesLo= [];

% Estimate the cost of each model
cost= zeros(1,np+1);
for j= 1:np+1
    img.params_mapping.params= modeles(:,j);
    cost(j)= objectiveFunction(img,data);
end
dist= 1;
% Proceed to the downhill simplex regression while reducing the temperature
while temp>=tempfinale && dist >= 1e-5
    for i= 1:nMetro
        [costlo,ilo]= min(cost);
        [costhi,ihi]= max(cost);
        [costnhi]= max(setdiff(cost,costhi));
        bary= sum(modeles,2);
        % Reflexion
        [modtryReflection,costtryReflection]= deformation(modeles,bary,np,ihi,-1,img,data);
        boolean1= annealingProbability(costtryReflection,costlo,temp);
        boolean2= annealingProbability(costtryReflection,costnhi,temp);
        if  boolean1 %costtryReflection <= costlo
             % Dilatation
            [modtryExpansion,costtryExpansion]= deformation(modeles,bary,np,ihi,-2,img,data);
            boolean= annealingProbability(costtryExpansion,costtryReflection,temp);
            if boolean %costtryExpansion <= costtryReflection
                modeles(:,ihi)= modtryExpansion;
                cost(ihi)= costtryExpansion;
            else
                modeles(:,ihi)= modtryReflection;
                cost(ihi)= costtryReflection; 
            end
        elseif boolean2 % costtryReflection<=costnhi
             % Outward  contraction
            [modtryOutContraction,costtryOutContraction]= deformation(modeles,bary,np,ihi,-0.7,img,data);
            boolean= annealingProbability(costtryOutContraction,costtryReflection,temp);
            if boolean % costtryOutContraction<=costtryReflection
                modeles(:,ihi)= modtryOutContraction;
                cost(ihi)= costtryOutContraction;
            else
                modeles(:,ihi)= modtryReflection;
                cost(ihi)= costtryReflection;     
            end     
        else 
            % Inward contraction
            [modtryContraction,costtryContraction]= deformation(modeles,bary,np,ihi,0.7,img,data);
            boolean= annealingProbability(costtryContraction,costhi,temp);
            if boolean % costtryContraction <= costhi
                modeles(:,ihi)= modtryContraction;
                cost(ihi)= costtryContraction;
             else   % Nothing better -> shrinkage
                modeles= (modeles+repmat(modeles(:,ilo),1,np+1))/2;
                for j= 1:np+1
                    img.params_mapping.params= modeles(:,j);
                    cost(j)= objectiveFunction(img,data);
                end
            end
        end
        costhi= max(cost);
        modelesT= [modelesT  modeles];
        costhhi= [costhhi ; costhi]; 
        [costlo,ilo]= min(cost);
        modelesLo= [modelesLo modeles(:,ilo)];
    end 
%      bary1= sum(modelesT(:,end-np:end),1)/(np+1);
%      bary2= sum(modelesT(:,end-2*(np+1)+1:end-(np+1)),1)/(np+1);
%      
%      dist= sqrt(sum((bary1-bary2).^2));
%      D(k)= dist;
     temperature(k)= temp;

    temp= temp*cooldown;
    niter= niter-1;
    k=k+1;
%     disp(['Remaining iterations= ' num2str(niter) ' - Temperature= ' num2str(temp)])

    
%     disp(['Remaining iterations= ' num2str(niter) ' - Temperature= ' num2str(temp) ' - Distance from previous iteration= ' num2str(dist)])
   
    if k>2 && D(k)>D(k-1); break; end
end
% modelesT1= reshape(modelesT(:,1),3,[]);
% modelesT2= reshape(modelesT(:,2),3,[]);


% figure; loglog(temperature,D); axis tight
% xlabel('Temperature','fontsize',20,'fontname','Times')
% ylabel('Models contraction','fontsize',20,'fontname','Times')
% set(gca,'fontsize',15,'fontname','Times','xdir','reverse')
% drawnow

[costlo,ilo]= min(cost);
img.params_mapping.params= modeles(:,ilo);
img= feval(mapping_function,img);


modeles(1:end-2,ilo);
exp(modeles(end-1:end,ilo));

v= 1:size(modelesLo,2);

figure; plot(v,exp(modelesLo),'linewidth',2);
xlabel('Iteration number','fontsize',20,'fontname','Times')
ylabel('Inversion parameter','fontsize',20,'fontname','Times');
set(gca,'fontsize',15,'fontname','Times');
axis tight;



figure; semilogy(reshape(costhhi,[],1),'k','linewidth',2)
xlabel('Iteration number','fontsize',20,'fontname','Times')
ylabel('Cost function','fontsize',20,'fontname','Times')
set(gca,'fontsize',15,'fontname','Times'); drawnow


end

function [boolean,Pr]= annealingProbability(costtry,costhi,temp)
E= (costtry-costhi)/temp;
Ptry= exp(-E);
Pr= min([1 Ptry]);
boolean= 0;
if Pr > rand ; boolean= 1; end
end


function cost= objectiveFunction(img,data)
vsim=  fwd_solve(img);
residuals= img.parameters.normalisation*(vsim.meas-data);
cost= sqrt(sum(residuals.^2));
end

function [modtry,costtry]= deformation(modeles,bary,np,ihi,fac,img,data)

fac1= (1-fac)/np;  % division by np since it is the weight of the barycenter
fac2= fac1-fac;      % no division as weight of a single model
% modification of the worst model
modtry= bary*fac1-modeles(:,ihi)*fac2;

imgtry= img;
imgtry.params_mapping.params= modtry;
costtry= objectiveFunction(imgtry,data);    
end



function do_unit_test
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;410,10,100) -maxh=20.0;\n'];
e0 = linspace(0,310,64)';
elec_pos = [e0,0*e0,0*e0,1+0*e0,0*e0,0*e0];
elec_shape= [0.1,0.1,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
fmdl.stimulation= stim_pattern_geophys( 64, 'Wenner', {'spacings', 1:32} );

img= mk_image(fmdl,1/20);
fm_pts= interp_mesh(fmdl);
x_bary= fm_pts(:,1); z_bary= fm_pts(:,2);

z_params= (min(fmdl.nodes(:,2)):max(fmdl.nodes(:,2)))';
a = 0.36; b = 130;
x_params= a*z_params+b;
xlim=interp1(z_params,x_params,z_bary);
img.elem_data(x_bary>xlim)= 1/120;

dd  = fwd_solve(img);
sig= sqrt(norm(dd.meas)); m= size(dd.meas,1);
noise= .05;
ddn= dd;
ddn.meas = dd.meas + noise*sig*randn(m,1);

img1= mk_image(fmdl,1);
vh1= fwd_solve(img1);
normalisation= 1./vh1.meas;
I= speye(length(normalisation));
I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

a = 0.3; b = 150;
res_params= log([10 100]');

z_params= (min(fmdl.nodes(:,2)):max(fmdl.nodes(:,2)))';
x_params= a*z_params+b;

imdl = eidors_obj('inv_model','testNoisy');
imdl.fwd_model= fmdl;
imdl.params_mapping.params= [a; b ; res_params];
imdl.params_mapping.perturb= [0.1; 50; 2 ; 2];

imdl.params_mapping.function = @border_mapping;
imdl.params_mapping.data.x_bary = x_bary;
imdl.params_mapping.data.z_bary = z_bary;
imdl.params_mapping.data.res_params = res_params;
imdl.params_mapping.data.x_params = x_params;
imdl.params_mapping.data.z_params = z_params;
imdl.params_mapping.data.a = a;
imdl.params_mapping.data.b = b;
imdl.reconst_type= 'absolute';
imdl.solve = @inv_solve_abs_annealingSimplex_params;
imdl.normalize_measurements= 1;
imdl.parameters.normalisation= I;

imdl.parameters.temp= 1000;
imdl.parameters.tempfinale= 1;
imdl.parameters.cooldown= 0.97;% 0.95;
imdl.parameters.nMetro= 1;
imdl.jacobian_bkgnd.value = 1;
% imgr= inv_solve(imdl, dd);
imgr= inv_solve(imdl, ddn);
img = mk_image( imdl );
img.elem_data= imgr.elem_data;
vAS= fwd_solve(img); vAS = vAS.meas;

figure; hist(I*(dd.meas-vAS),50)
show_pseudosection( fmdl, I*dd.meas, 'HorizontalDownward')
show_pseudosection( fmdl, I*vAS, 'HorizontalDownward')
show_pseudosection( fmdl, (vAS-dd.meas)./dd.meas*100,'HorizontalDownward')

end


function img = border_mapping(img)
%% Function to be called to perform the mapping in the forward problem
z= img.params_mapping.data.z_params;
res= img.params_mapping.params(end-1:end);

a= img.params_mapping.params(1);
b= img.params_mapping.params(2);
x= z*a+b; 

xi= img.params_mapping.data.x_bary;
zi= img.params_mapping.data.z_bary;
xlim=interp1(z,x,zi);

vi= zeros(size(img.fwd_model.elems,1),1) + res(1);
vi(xi>xlim)= res(2);

img.elem_data= exp(-vi);

img.params_mapping.data.x_params = x;
img.params_mapping.params= [a ; b ;  res];

end

