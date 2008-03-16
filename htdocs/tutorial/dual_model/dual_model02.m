% Simulate data $Id: dual_model02.m,v 1.1 2008-03-16 00:12:52 aadler Exp $

imdl=mk_common_model('a2c0',16);

for model = 1:2
   if model==1
% Model 1: coarse==fine. each elem has a parameter
      params= 1:64; 
   else
% Model 2: outer two layers have only one parameter, 
      params= 1:64;
      params(params>36)= 37;
   end

   subplot(2,2, model)
   show_fem(imdl.fwd_model);

% Show parameter numbers
   elems= imdl.fwd_model.elems;
   nodes= imdl.fwd_model.nodes;
   e= size(elems,1);
   numeros= reshape(sprintf('%2d',params),2,e)';
   xc=mean(reshape(nodes(elems,1),e,3),2);
   yc=mean(reshape(nodes(elems,2),e,3),2);
   text(xc,yc,numeros,'FontSize',8, ...
            'HorizontalAlignment','center');

end

print -r125 -dpng dual_model02a.png;
