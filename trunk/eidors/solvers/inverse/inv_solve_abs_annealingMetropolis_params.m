function [img]= inv_solve_abs_annealingMetropolis_params(inv_model, data)
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
% $Id: CG_lr_solve.m 2596 2012-06-06 14:18:06Z nlesparre $


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
    nMetro= 300;
end

if ~isfield(inv_model.parameters,'normalisation')
    img.parameters.normalisation= 1;
end

temp= tempInit; k= 1;
while temp>=tempfinale
    temp= temp*cooldown;
    k=k+1;
end

disp(['Number of iteration ' num2str(k)])

% Estimate the number of parameters to adjust
n_parameters= size(img.params_mapping.params,1);

% Generates np+1 models with a Gaussian law with a mean equals to
% inv_model.params_mapping.params and a standard deviation of
% inv_model.params_mapping.perturb
if isfield(inv_model.params_mapping,'inital_model')
    model= inv_model.params_mapping.inital_model;
else
    range= inv_model.params_mapping.range;
    model= 10.^((rand(np,1).*(range(:,2)-range(:,1))+range(:,1)));    
end

% Estimate the cost of the initial model
img.params_mapping.params= model;
cost= objectiveFunction(img,data);


% Proceed to the Metropolis regression while reducing the temperature
%% Run over the Simulated annealing loops

temperature= t_start;
n_temp= 1;
while temperature >= t_end % Loop over temperature
    for n= 1:nMetro % Metropolis loop
        model_try= model;
        idx_parametr= randi(n_parameters,1,1);
        model_try(idx_parametr)= 10.^(rand(1,1)*(range(idx_parametr,2)-range(idx_parametr,1))+range(idx_parametr,1));
        img.params_mapping.params= model_try;
        cost_try= objectiveFunction(img,data);
        if rand <= (cost_try/cost)^(1/temperature) % Metropolis test for acceptance of trial model
            likelihood= (cost_try/cost)^(1/temperature);
            cost= cost_try;
            model= model_try;
        end
    end % End of Metropolis loop
    resRec(n_temp)= cost;
    res_parRec(n_temp,:)= model';
    tempRec(n_temp)= temperature;
    likelihoodRec(n_temp)= likelihood;
        
    n_temp= n_temp+1;
    temperature = temperature*alpha;	% Decrease the temperature
end % End of temperature loop

screenSize = get(0,'ScreenSize');
h = figure;
set(h,'Position',[0 0 screenSize(3)/2 screenSize(4)]);
subplot(4,1,1)
plot(tempRec,res_parRec(:,1))
set(gca,'fontsize',16,'fontname','Times','xDir','reverse','xScale','log','yScale','log')
subplot(4,1,2)
plot(tempRec,res_parRec(:,2))
set(gca,'fontsize',16,'fontname','Times','xDir','reverse','xScale','log','yScale','log')
subplot(4,1,3)
plot(tempRec,resRec)
set(gca,'fontsize',16,'fontname','Times','xDir','reverse','xScale','log','yScale','log')
subplot(4,1,4)
plot(tempRec,likelihoodRec)
set(gca,'fontsize',16,'fontname','Times','xDir','reverse','xScale','log')

end

function cost= objectiveFunction(img,data)
vsim=  fwd_solve(img);
residuals= img.parameters.normalisation*(vsim.meas-data);
cost= sqrt(sum(residuals.^2));
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

