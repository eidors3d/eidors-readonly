% load some lung data
load montreal_data_1995;
vi = zc_resp(:,2); 
vh = zc_resp(:,19); 

% Lung forward model - with electrode #1 on back
imdl_lung= mk_common_model('b2t2',16);
fmdl_lung= imdl_lung.fwd_model;
fmdl_lung.electrode= fmdl_lung.electrode([9:16,1:8]); 

% Put lung models into the inv_model for elems
imdl_e.fwd_model= fmdl_lung;

% Put lung models into the inv_model for nodes
fmdl_lung.jacobian_elem2nodes.fwd_model = fmdl_lung;
fmdl_lung.jacobian = @jacobian_elem2nodes;
imdl_n.fwd_model= fmdl_lung;
imdl_n.RtR_prior_elem2nodes.fwd_model = fmdl_lung;

img1= inv_solve(imdl_e,vh,vi);
subplot(121)
show_slices(img1);
img2= inv_solve(imdl_n,vh,vi);
subplot(122)
show_slices(img2);


%print -r125 -dpng nodal_jacobian_solver.png;

