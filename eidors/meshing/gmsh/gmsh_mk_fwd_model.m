function [fwd_mdl, mat_indices]= ...
             gmsh_mk_fwd_model( vol_filename, name, centres, z_contact)
% GMSH_MK_FWD_MODEL: create a fwd_model object from a Gmsh file
% [fwd_mdl, mat_indices]= ...
%      gmsh_mk_fwd_model( vol_filename, centres, ...
%                       name, stim_pattern, z_contact)
%
%  vol_filename:      filename output from Gmsh
%  name:              name for object (if [] use vol_filename)
%  centres:           matrix of N x [x,y,z] electrode centres
%                     centres can also be a Nx1 cell matrix of
%                     functions which are 1 inside the electrode and 0 outside
%  stim_pattern:      a stimulation pattern structure
%                     empty ([]) if stim_pattern is not available
%  z_contact:         vector or scalar electrode contact impedance
%
%  fwd_mdl:           eidors format fwd_model
%  mat_indices:       cell array of material indices from eidors 

% Gmsh mesher for EIRODS was based on Netgen interface.
% (C) 2009 Bartosz Sawicki. License: GPL version 2 or version 3

if isempty(name); 
   name = ['fwd_mdl based on ', vol_filename];
end

if nargin<5
   z_contact=0.01; % singular if z_contact=0
end

stim_pattern=[];
% Model Geometry
[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh(vol_filename);
fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                             stim_pattern, centres, z_contact);
mat_indices= mk_mat_indices( mat_ind);
if isempty(srf)
   fwd_mdl.boundary = find_boundary(fwd_mdl);
end

fwd_mdl.mat_idx = mat_indices;

% build fwd_model structure
function fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
    stim_pattern, centres, z_contact)
    mdl.nodes    = vtx;
    mdl.elems    = simp;
    mdl.boundary = srf;
    mdl.gnd_node = 1;
    mdl.name = name;
    
    % Model Stimulation
    if ~isempty(stim_pattern)
        mdl.stimulation= stim_pattern;
    end
    
    nelec= size(centres,1);
    % Electrodes
%    [elec,sels,electrodes] = ng_tank_find_elec(srf,vtx,bc,centres);
%    if size(elec,1) ~= nelec
%        error('Failed to find all the electrodes')
%    end
    

    % set the z_contact
    z_contact= z_contact.*ones(nelec,1);
    for i=1:nelec
        electrodes(i).nodes(1) = i; 
        electrodes(i).z_contact= z_contact(i);
    end
    if nelec >0
    mdl.electrode =     electrodes;
    end
    mdl.solve=          'eidors_default';
    mdl.jacobian=       'eidors_default';
    mdl.system_mat=     'eidors_default';

    fwd_mdl= eidors_obj('fwd_model', mdl);

% Output cell array of indices into each material type
%   array order is sorted by length of material type
function mat_indices= mk_mat_indices( mat_ind);
    % find length of mat_indices 
    % test example: mat_ind=[10 12 14 14 12 12 14 12];
    sort_mi= sort(mat_ind(:));
    find_mi= find( diff([-1e8;sort_mi]) );
    len_mi = diff([find_mi;length(sort_mi)+1]);
    [jnk,idxs]= sort(-len_mi); %reverse sort
    l_idxs= length(idxs);
    mat_indices= cell(1, l_idxs);
    for i= 1:l_idxs;
        mat_idx_i= sort_mi(find_mi(idxs(i)));
        mat_indices{i}= find(mat_ind == mat_idx_i);
    end

