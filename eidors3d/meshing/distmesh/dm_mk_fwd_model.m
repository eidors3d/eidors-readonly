function [fwd_mdl]= dm_mk_fwd_model( fd, fh, h0, bbox, ...
                               stim_pattern, z_contact, name)
% DM_MK_FWD_MODEL: create a fwd_model object using distmesh
% fwd_mdl= dm_mk_fwd_model( fd, fh, h0, bbox,...
%                          elec_nodes, stim_pattern, z_contact, name);
%
%  FD:        Distance function - 1 for space inside area, 0 outside
%  FH:        Scaled edge length function - decrease in refined areas
%  H0:        Initial edge length
%  BBOX:      Bounding box [xmin,ymin, {zmin}; xmax,ymax,{zmax}]
%
%  centres:           vector of electrode centres from 'create_tank_mesh_ng' 
%  stim_pattern:      a stimulation pattern structure
%                     empty ([]) if stim_pattern is not available
%  z_contact:         vector or scalar electrode contact impedance
%  name:              name for eidors object
%
%  fwd_mdl:           eidors format fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_mk_fwd_model.m,v 1.2 2008-03-10 14:41:28 aadler Exp $

if nargin <8
   name = 'MDL from dm_mk_fwd_model';
end

fix_node= [];
keyboard
[vtx,simp] = distmeshnd(fd,fh,h0,bbox,fix_node);
srf= find_boundary(simp);

fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                       stim_pattern, centres, z_contact)

% build fwd_model structure
function fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                       stim_pattern, centres, z_contact)
mdl.nodes    = vtx;
mdl.elems    = simp;
mdl.boundary = srf;
mdl.gnd_node=           1;
mdl.np_fwd_solve.perm_sym =     '{n}';
mdl.name = name;

% Model Stimulation
if ~isempty(stim_pattern)
   mdl.stimulation= stim_pattern;
end

% Electrodes

% set the z_contact
z_contact= z_contact.*ones(nelec,1);
for i=1:nelec
   electrodes(i).z_contact= z_contact(i);
end

mdl.electrode =     electrodes;
mdl.solve=          'np_fwd_solve';
mdl.jacobian=       'np_calc_jacobian';
mdl.system_mat=     'np_calc_system_mat';

fwd_mdl= eidors_obj('fwd_model', mdl);
