function [fwd_mdl, mat_indices]= ...
             ng_mk_fwd_model( ng_vol_filename, centres, ...
                              name, stim_pattern, z_contact, postprocmesh)
% NG_MK_FWD_MODEL: create a fwd_model object from a netgen vol file
% [fwd_mdl, mat_indices]= ...
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
%  mat_indices:       cell array of material indices from eidors 
 
% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isempty(name); 
   name = ['MDL from', ng_vol_filename];
end

if nargin<5
   z_contact=0.01; % singular if z_contact=0
end

% Model Geometry
[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(ng_vol_filename);
if nargin>=6
    N_elec = max(size(centres));
    [srf,vtx,fc,bc,simp,edg,mat_ind] = feval(postprocmesh,...
        srf,vtx,fc,bc,simp,edg,mat_ind, N_elec);
end
fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                             stim_pattern, centres, z_contact,fc);
mat_indices= mk_mat_indices( mat_ind);

% build fwd_model structure
function fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                       stim_pattern, centres, z_contact,fc)
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

mdl.solve=      @fwd_solve_1st_order;
mdl.jacobian=   @jacobian_adjoint;
mdl.system_mat= @system_mat_1st_order;

fwd_mdl= eidors_obj('fwd_model', mdl);

% Output cell array of indices into each material type
%   array order is the order of the specified material,
%   except that the largest material (background) is placed
%   first.
function mat_idx= mk_mat_indices( mat_ind);
  % find length of mat_indices 
  % test example: mat_ind=[10 12 14 14 12 12 14 12];

  mat_indices = unique( mat_ind );
  for i= 1:length(mat_indices);
     mat_idx{i}= find(mat_ind == mat_indices(i));
     mat_idx_l(i) = length( mat_idx{i} );
  end
  % put the largest material first
  [jnk, max_l] = max(mat_idx_l);
  new_idx = 1:length(mat_indices); new_idx(max_l) = [];
  ver= eidors_obj('interpreter_version');
  if ~ver.isoctave && ver.ver < 7
% STUPID MATLAB WORK AROUND
     clear jnk;ii=1;for i=[max_l, new_idx];
         jnk{ii} = mat_idx{i}; ii=ii+1;
     end
     mat_idx = jnk;
  else
     [mat_idx{:}] = mat_idx{[max_l, new_idx]};
  end


function gnd_node=    find_centre_node(vtx);
  %distance from zero
  d = sum( vtx.^2, 2);
  [jnk,gnd_node] = min(d);
  gnd_node= gnd_node(1);
