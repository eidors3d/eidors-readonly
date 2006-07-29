function fwd_mdl= ng_mk_fwd_model( ng_vol_filename, centres, ...
                                   name, stim_pattern, z_contact)
% NG_MK_FWD_MODEL: create a fwd_model object from a netgen vol file
% fwd_mdl= ng_mk_fwd_model( ng_vol_filename, name, ...
%                           stim_pattern, z_contact)
%
%  ng_vol_filename:   filename output from netgen
%  name:              name for object (if [] use ng_vol_filename)
%  centres:           vector of electrode centres from 'create_tank_mesh_ng' 
%  stim_pattern:      a stimulation pattern structure
%                     empty ([]) if stim_pattern is not available
%  z_contact:         vector or scalar electrode contact impedance
%
% (C) 2006 Andy Adler. Licenced under GPL V2

if nargin<5
   z_contact=0;
end

mdl.name = name;
if isempty(name); 
   mdl.name = ['MDL from', ng_vol_filename];
end

% Model Geometry
[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(ng_vol_filename);
mdl.nodes    = vtx;
mdl.elems    = simp;
mdl.boundary = srf;
mdl.gnd_node=           1;
mdl.misc.perm_sym =     '{n}';

% Model Stimulation
if ~isempty(stim_pattern)
   mdl.stimulation= stim_pattern;
end

nelec= size(centres,1);
% Electrodes
[elec,sels] = ng_tank_find_elec(srf,vtx,bc,centres);
if size(elec,1) ~= nelec
   error('Failed to find all the electrodes')
end

z_contact= z_contact.*ones(nelec,1);
for i=1:nelec
   electrodes(i).z_contact= z_contact(i);
   electrodes(i).nodes=     unique( elec(i,:) );
end

mdl.electrode =         electrodes;
mdl.solve=          'np_fwd_solve';
mdl.jacobian=       'np_calc_jacobian';
mdl.system_mat=     'np_calc_system_mat';

fwd_mdl= eidors_obj('fwd_model', mdl);
