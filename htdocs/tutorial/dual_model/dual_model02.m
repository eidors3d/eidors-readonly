% Simulate data $Id$

% create base model
mdl_base=mk_common_model('a2c0',16);
mdl_base.RtR_prior = @noser_image_prior;
mdl_base.hyperparameter.value = 3e-2;

elems= mdl_base.fwd_model.elems;
nodes= mdl_base.fwd_model.nodes;
e= size(elems,1);


for model = 1:3
   if model==1
% Model 1: coarse==fine. each elem has a parameter
      params= 1:e; 
   elseif model==2
% Model 2: coarse model, inner circle has one parameter
      params= [1,1,1,1, 2:e-3];
   elseif model==3
% Model 3: coarse model, top left slice has one parameter
      params= 1:e;
      params([4,8,15:16,23:24,34:36])= 0;
      [jnk1,jnk2,params]= unique(params);
   end

% Create inverse_model
   imdl(model)= mdl_base;
   imdl(model).fwd_model.coarse2fine = sparse(1:e,params,1,e,max(params));

   subplot(2,3, model)
   show_fem(imdl(model).fwd_model);

% Show parameter numbers
   numeros= reshape(sprintf('%2d',params),2,e)';
   xc=mean(reshape(nodes(elems,1),e,3),2);
   yc=mean(reshape(nodes(elems,2),e,3),2);
   text(xc,yc,numeros,'FontSize',8, ...
            'HorizontalAlignment','center');

end

print -r125 -dpng dual_model02a.png;
