% Make square mesh $Id$

% Create square mesh model
imdl= mk_common_model('c2s',16);
s_mdl= rmfield(imdl.fwd_model,{'electrode','stimulation'});

% assign one parameter to each square
e= size(s_mdl.elems,1);
params= ceil(( 1:e )/2);
s_mdl.coarse2fine = sparse(1:e,params,1,e,max(params));

show_fem(s_mdl)

% Show parameter numbers
   numeros= reshape(sprintf('%3d',params),3,e)';
   xc=mean(reshape(s_mdl.nodes(s_mdl.elems,1),e,3),2);
   yc=mean(reshape(s_mdl.nodes(s_mdl.elems,2),e,3),2);
   text(xc,yc,numeros,'FontSize',7, ...
            'HorizontalAlignment','center');

print -r100 -dpng square_mesh01a.png;
