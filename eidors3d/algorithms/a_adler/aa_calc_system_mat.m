function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_calc_system_mat.m,v 1.17 2008-05-11 00:27:13 aadler Exp $

FC= aa_system_mat_fields( fwd_model);

pp= aa_fwd_parameters( fwd_model);
elem_sigma = kron( [img.elem_data(:);ones(pp.n_elec,1)], ...
                    ones(pp.n_dims,1) );
l_elem_s= length( elem_sigma);
elem_sigma = spdiags( elem_sigma, 0, l_elem_s, l_elem_s);

s_mat.E= FC' * elem_sigma * FC;

