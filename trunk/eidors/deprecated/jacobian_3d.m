function [J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,df,tol,perm_sym);
%function [J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,df,tol,perm_sym);
%
%This function calculates the Jacobian (sensitivity) matrix of the system at
%mat_ref. 
%
%
%
%I        = The currents used
%elec     = the electrodes matrix
%vtx      = The vertices matrix
%simp     = The simplices matrix
%gnd_ind  = The ground index (node)
%mat_ref  = The reference conductivity vector
%zc       = The electrode contact impedance vector
%IntGrad  = The integrals of the gradients
%v_f      = The measurement fields
%df       = Measurements per current pattern as used in v_f 
%tol      = Tolerance 
%J        = The Jacobian (sensitivity) matrix with respect to conductivity

warning('EIDORS:deprecated','JACOBIAN_3D is deprecated as of 07-Jun-2012. ');

[vr,vc] = size(vtx);
[sr,sc] = size(simp);

el_no = size(elec,1);

if sum(df)~= size(v_f,2);
   error('Mismatched data input');
end

[E,D,Ela,pp] = fem_master_full(vtx,simp,mat_ref,gnd_ind,elec,zc,perm_sym);

[V] = forward_solver(E,I,tol,pp);

[J] = jacobian_3d_fields(V,Ela,D,elec,vtx,simp,mat_ref,v_f,df);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides and W.R.B. Lionheart 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
