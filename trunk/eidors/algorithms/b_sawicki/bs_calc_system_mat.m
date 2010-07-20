function s_mat= bs_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% Modified by Bartosz Sawicki
% $Id: bs_calc_system_mat.m 1535 2008-07-26 15:36:27Z aadler $

FC= aa_system_mat_fields( fwd_model);
lFC= size(FC,1);

pp= aa_fwd_parameters( fwd_model);

if( size(img.elem_data,1) == pp.n_elem )
    %For non-dual models
    elem_data = img.elem_data; 
else
    %For dual models
    if isfield(fwd_model, 'coarse2fine')
        elem_data = fwd_model.coarse2fine * img.elem_data;
    end  
    %If background image is known
    if isfield(fwd_model, 'background')
        elem_data = elem_data + fwd_model.background;
    end
end
  
elem_sigma = kron( elem_data(:), ones(pp.n_dims,1) );

ES= ones(lFC,1);
ES(1:length(elem_sigma))= elem_sigma;
ES= spdiags(ES,0,lFC,lFC);

s_mat.E= FC' * ES * FC;

