function s_mat= system_mat_1st_order( fwd_model, img)
% SYSTEM_MAT_1ST_ORDER: SS= system_mat_1st_order( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

FC= system_mat_fields( fwd_model);
lFC= size(FC,1);

elem_sigma = kron( img.elem_data(:), ones(elem_dim(fwd_model),1) );

ES= ones(lFC,1);
ES(1:length(elem_sigma))= elem_sigma;
ES= spdiags(ES,0,lFC,lFC);

s_mat.E= FC' * ES * FC;

