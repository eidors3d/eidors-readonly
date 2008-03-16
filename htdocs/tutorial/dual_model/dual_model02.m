% Simulate data $Id: dual_model02.m,v 1.3 2008-03-16 10:35:23 aadler Exp $

% create base model
mdl_base=mk_common_model('a2c0',16);
mdl_base.RtR_prior = @noser_image_prior;
mdl_base.hyperparameter.value = 3e-2;

elems= mdl_base.fwd_model.elems;
nodes= mdl_base.fwd_model.nodes;
e= size(elems,1);


for model = 1:2
   if model==1
% Model 1: coarse==fine. each elem has a parameter
      params= 1:e; 
   else
% Model 2: outer two layers have only one parameter, 
      params= 1:e;
      params(params>36)= 37;
   end

% Create inverse_model
   imdl(model)= mdl_base;
   imdl(model).fwd_model.coarse2fine = sparse(1:e,params,1,e,max(params));

   subplot(2,2, model)
   show_fem(imdl(model).fwd_model);

% Show parameter numbers
   numeros= reshape(sprintf('%2d',params),2,e)';
   xc=mean(reshape(nodes(elems,1),e,3),2);
   yc=mean(reshape(nodes(elems,2),e,3),2);
   text(xc,yc,numeros,'FontSize',8, ...
            'HorizontalAlignment','center');

end

print -r125 -dpng dual_model02a.png;
