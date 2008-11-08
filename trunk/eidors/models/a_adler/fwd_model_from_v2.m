function fmdl = fwd_model_from_v2( elec, vtx, bdy, simp, zc, gnd_ind, name)
% CREATE EIDORS v3 fwd_model from the v2 variables
% fmdl = fwd_model_from_v2( elec, vtx, simp, zc, gnd_ind)
% Input paramers
%  elec: [N_elec x 3*Srf] for each electrode
%  vtx:  [N_nodes x 3] Vertices
%  bdy:  [N_bdy x 2]   Boundary elements (specify [] if not available)
%  simp:  [N_elems x 4] Tetrahedral elements
%  zc;   [N_elec x 1]  Electrode contact impedances
%  gnd_ind: [1x1] ground Vertex
%  name:  fwd_model name (optional)
% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isempty(bdy);
   bdy= find_boundary( simp );
end

fmdl.type =     'fwd_model';
fmdl.name =     'fwd_model_from_v2';
fmdl.nodes=     vtx;
fmdl.elems=     simp;
fmdl.boundary=  bdy;
fmdl.electrode= get_model_elecs(elec,zc);
fmdl.solve=     'np_fwd_solve';
fmdl.jacobian=  'np_calc_jacobian';
fmdl.system_mat='np_calc_system_mat';
fmdl.gnd_node=  gnd_ind;
fmdl.np_fwd_solve.perm_sym = '{n}'; % default, can change


function [electrodes] = get_model_elecs(elec,zc);
%elec : The electrodes matrix. 
%np_pl : Number of electrode planes (in planar arrangements)
%zc : Contact impedances of the electrodes

for i=1:length(zc)
    electrodes(i).z_contact= zc(i);
    electrodes(i).nodes=     unique( elec(i,:) );
end

perm_sym='{n}';
