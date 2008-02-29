%$Id: nodal_jacobian01.m,v 1.2 2008-02-29 22:50:43 aadler Exp $
imdl_e= mk_common_model('b2t2',16);
imdl_e.hyperparameter.value= 0.05;
imdl_e.RtR_prior = @laplace_image_prior;

fmdl= imdl_e.fwd_model;
imdl_n= imdl_e;
imdl_n.fwd_model.jacobian = @jacobian_elem2nodes;
imdl_n.fwd_model.jacobian_elem2nodes.fwd_model = fmdl;
imdl_n.RtR_prior = @RtR_prior_elem2nodes;
imdl_n.RtR_prior_elem2nodes.RtR_prior = @laplace_image_prior;
imdl_n.RtR_prior_elem2nodes.fwd_model = fmdl;
imdl_n.reconst_to = 'nodes';

img_e= calc_jacobian_bkgnd( imdl_e );
img_n= img_e;
img_n.fwd_model= imdl_n.fwd_model;
img_n= rmfield(img_n,'elem_data');
img_n.node_data= mapper_nodes_elems(fmdl)*img_e.elem_data;


J_e= calc_jacobian( img_e);
J_n= calc_jacobian( img_n);

% load some lung data
load montreal_data_1995;
vi = zc_resp(:,2); 
vh = zc_resp(:,19); 

img1= inv_solve(imdl_e,vh,vi);
show_slices(img1);
img2= inv_solve(imdl_n,vh,vi);
NtoE = mapper_nodes_elems( fmdl );


%print -r125 -dpng nodal_jacobian_solver.png;
