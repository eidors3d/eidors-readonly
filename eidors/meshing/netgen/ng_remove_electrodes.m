function [srf,vtx,fc,bc,simp,edg,mat_ind] = ng_remove_electrodes...
    (srf,vtx,fc,bc,simp,edg,mat_ind, N_elec)
% NG_REMOVE_ELECTRODES: cleans up matrices read from a *.vol file
% [srf,vtx,fc,bc,simp,edg,mat_ind]= ng_remove_electrodes...
%     (srf,vtx,fc,bc,simp,edg,mat_ind, N_elec)
%
% Used to clean up external objects used to force electrode meshing in
% ng_mk_extruded_model.
%
% (C) Bartlomiej Grychtol, 2010. Licenced under GPL v2 or v3
% $Id$

% total objects:
N_obj = max(mat_ind);

% electrode simps:
e_simp_ind = mat_ind > (N_obj - N_elec);

in = unique(simp(~e_simp_ind,:));
out = unique(simp(e_simp_ind,:));
boundary = intersect(in,out);
out = setdiff(out,boundary);

ext_srf_ind = ismember(srf,out);
ext_srf_ind = ext_srf_ind(:,1) | ext_srf_ind(:,2) | ext_srf_ind(:,3);

srf(ext_srf_ind,:) = [];
bc(ext_srf_ind,:) = [];
fc(ext_srf_ind,:) = [];
simp = simp(~e_simp_ind,:);
mat_ind = mat_ind(~e_simp_ind);

% fix bc:
n_unique = numel(unique(bc));
missing = setdiff(1:n_unique, unique(bc));
spare = setdiff(unique(bc), 1:n_unique); 
for i = 1:length(missing)
    bc( bc==spare(i) ) = missing(i);
end

% fic vtx:
v = 1:size(vtx,1);
unused_v = setdiff(v, union(unique(simp),unique(srf))); 
v(unused_v) = [];
for i = 1: length(vtx)
%     simp_ind = find(simp == i);
%     srf_ind = find( srf == i);
    new_v_ind = find(v == i);
    simp( simp == i ) = new_v_ind; 
    srf( srf  == i ) = new_v_ind;
end
vtx(unused_v,:) = [];


