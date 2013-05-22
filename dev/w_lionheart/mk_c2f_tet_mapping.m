function mapping = mk_c2f_tet_mapping( f_mdl, c_mdl );
% MK_C2F_TET_MAPPING: create a mapping matrix from coarse to fine FEM
% c2f= mk_c2f_tet_mapping( f_mdl, c_mdl );
%  
% Parameters:
%    c_mdl is coarse fwd_model
%    f_mdl is fine fwd_model
%
% C2F_ij is the fraction if f_mdl element i which is
%   contained in c_mdl element j. This is used to map
%   from data on the reconstruction model (c_mdl) to
%   the forward model f_mdl as 
%      elem_data_fine = Mapping*elem_data_coase

if isstr(f_mdl) && strcmp(f_mdl, 'UNIT_TEST'); do_unit_test; return; end
