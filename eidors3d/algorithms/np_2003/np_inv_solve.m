function img= np_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% $Id: np_inv_solve.m,v 1.4 2004-07-18 03:17:25 aadler Exp $

fwd_model= inv_model.fwd_model;
p= np_fwd_parameters( fwd_model );

% calc jacobian with homogeneous background
homg_img.type = 'image';
homg_img.elem_data= ones( p.n_elem ,1);
homg_img.fwd_model= fwd_model;

J = calc_jacobian( fwd_model, homg_img);

% Calculating a smoothing prior
[Reg] = iso_f_smooth(p.simp,p.vtx,3,1);

% Calculating a linear inverse solution

tfac= inv_model.hyperparameter;
dva= data1.meas - data2.meas;
sol = (J'*J +  tfac*Reg'*Reg)\J' * dva;


% create a data structure to return
img.name= 'solved by np_inv_solve';
img.elem_data = sol;
img.type = 'real conductivity differences';
img.fwd_model= fwd_model;
img.inv_model= inv_model;
