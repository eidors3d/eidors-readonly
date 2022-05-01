function [fwd_mdl]= ...
             ng_mk_fwd_model( ng_vol_filename, centres, ...
                              name, stim_pattern, z_contact, postprocmesh)
% NG_MK_FWD_MODEL: create a fwd_model object from a netgen vol file
% [fwd_mdl]= ...
%      ng_mk_fwd_model( ng_vol_filename, centres, ...
%                       name, stim_pattern, z_contact)
%
%  ng_vol_filename:   filename output from netgen
%  name:              name for object (if [] use ng_vol_filename)
%  centres:           matrix of N x [x,y,z] electrode centres
%                     centres can also be a Nx1 cell matrix of
%                     functions which are 1 inside the electrode and 0 outside
%  stim_pattern:      a stimulation pattern structure
%                     empty ([]) if stim_pattern is not available
%  z_contact:         vector or scalar electrode contact impedance
%
%  fwd_mdl:           eidors format fwd_model
%  fwd_mdl.mat_idx:   cell array of indices into netgen material properties
% (C) 2006 Andy Adler. (C) 2013 Alistair Boyle. License: GPL version 2 or version 3
% $Id$

if nargin==1 && strcmp(ng_vol_filename,'UNIT_TEST'); do_unit_test; return; end

if isempty(name);
   name = ['MDL from', ng_vol_filename];
end

if nargin<5
   z_contact=0.01; % singular if z_contact=0
end

% Model Geometry
[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(ng_vol_filename);
if isempty(vtx);
   error('EIDORS: ng_mk_fwd_model: Netgen meshing failed. Stopping');
end
if nargin>=6
    N_elec = size(centres,1);
    [srf,vtx,fc,bc,simp,edg,mat_ind] = feval(postprocmesh,...
        srf,vtx,fc,bc,simp,edg,mat_ind, N_elec);
end
fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                             stim_pattern, centres, z_contact,fc);

[fwd_mdl.mat_idx] = mk_mat_indices(mat_ind);

if ~isfield(fwd_mdl,'normalize_measurements')
   fwd_mdl.normalize_measurements = 0;
end

% build fwd_model structure
function fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                       stim_pattern, centres, z_contact,fc)
% maybe there were bugs from an unsorted boundary. Just in case

% get rid of duplicates, if there are any
[~,idx] = unique(sort(srf,2),'rows','stable');
if length(idx)~=size(srf,1)
   error('Duplicate elements on boundary. Netgen completed, but EIDORS can''t interpret the result. Try modifying the netgen parameters.')
end
%fc = fc(idx);
%bc = bc(idx);

mdl.nodes    = vtx;
mdl.elems    = simp;
mdl.boundary = srf;
mdl.boundary_numbers=fc;
mdl.gnd_node=    find_centre_node(vtx);
mdl.np_fwd_solve.perm_sym =     '{n}';
mdl.name = name;

% Model Stimulation
if ~isempty(stim_pattern)
   mdl.stimulation= stim_pattern;
end

nelec= size(centres,1);
if nelec>0
   % Electrodes
   [elec,sels,electrodes] = ng_tank_find_elec(srf,vtx,bc,centres);
   if size(elec,1) ~= nelec
      error('Failed to find all the electrodes')
   end

   % set the z_contact
   z_contact= z_contact.*ones(nelec,1);
   for i=1:nelec
      electrodes(i).z_contact= z_contact(i);
   end

   mdl.electrode =     electrodes;
end

mdl.solve=      'eidors_default';
mdl.jacobian=   'eidors_default';
mdl.system_mat= 'eidors_default';

fwd_mdl= eidors_obj('fwd_model', mdl);

% Output mat_idx cell array of indices into each material
% typei. Array order is the order of the specified material
% (netgen 'tlo' statements in the .geo file).
function [mat_idx] = mk_mat_indices( mat_ind);
  % find length of mat_indices
  % test example: mat_ind=[10 12 14 14 12 12 14 12];

  if isempty(mat_ind)
     mat_idx = [];
     return
  end
  mat_indices = unique( mat_ind );
  for i= 1:length(mat_indices);
     mat_idx{i}= find(mat_ind == mat_indices(i));
  end

function gnd_node=    find_centre_node(vtx);
  %distance from zero
  d = sum( vtx.^2, 2);
  [jnk,gnd_node] = min(d);
  gnd_node= gnd_node(1);

function do_unit_test
  % Test we can make one electrode
  fmdl = ng_mk_cyl_models(1,[],[]);
  unit_test_cmp('No elecs',num_elecs(fmdl),0);
  fmdl = ng_mk_cyl_models(1,[0,0.5;10,0.2],[0.1]);
  unit_test_cmp('Two elecs',num_elecs(fmdl),2);
