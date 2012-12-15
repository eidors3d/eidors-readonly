function [fwd_mdl, mat_indices]= ...
             gmsh_mk_fwd_model( vol_filename, name, eprefix,stim_pattern, z_contact)
% GMSH_MK_FWD_MODEL: create a fwd_model object from a Gmsh file
% [fwd_mdl, mat_indices]= ...
%      gmsh_mk_fwd_model( vol_filename, centres, ...
%                       name, stim_pattern, z_contact)
%
%  vol_filename:      filename output from Gmsh
%  name:              name for object (if [] use vol_filename)
%  eprefix:           prefix used for names of electrodes
%                     (if [] or omitted use 'electrode-')
%  stim_pattern:      a stimulation pattern structure
%                     empty ([]) if stim_pattern is not available
%  z_contact:         vector or scalar electrode contact impedance
%
%  fwd_mdl:           eidors format fwd_model
%  mat_indices:       cell array of material indices from eidors 

% Gmsh mesher for EIDORS was based on Netgen interface.
% (C) 2009 Bartosz Sawicki. License: GPL version 2 or version 3
% Modified by James Snyder & Bartlomiej Grychtol

if isempty(name); 
   name = ['fwd_mdl based on ', vol_filename];
end

if nargin < 4
    stim_pattern=[];
end

if nargin<3 || isempty(eprefix); 
   eprefix = 'electrode-';
end

if nargin<5
   z_contact=0.01; % singular if z_contact=0
end


% Model Geometry
[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(vol_filename);
fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                             stim_pattern, eprefix, z_contact, fc,phys_names);
mat_indices= mk_mat_indices( mat_ind);
if isempty(srf)
   fwd_mdl.boundary = find_boundary(fwd_mdl);
end

fwd_mdl.mat_idx = mat_indices;

% build fwd_model structure
function fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
    stim_pattern, eprefix, z_contact, fc, phys_names)
    mdl.nodes    = vtx;
    mdl.elems    = simp;
    mdl.boundary = srf;
    mdl.boundary_numbers=fc; 
    mdl.gnd_node = 1;
    mdl.name = name;
    
    % Model Stimulation
    if ~isempty(stim_pattern)
        mdl.stimulation= stim_pattern;
    end
    
    % Electrodes
    electrodes = find_elec(phys_names,eprefix,z_contact);
    if ~isempty(fields(electrodes));
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

% Assumes that electrodes are numbered starting at 1, with prefix provided
function electrodes = find_elec(phys_names,prefix,z_contact)
electrodes = struct();
phys_elecs = find(arrayfun(@(x)strncmp(x.name,prefix,length(prefix)),phys_names));
for i = 1:length(phys_elecs)
    cur_elec = arrayfun(@(x)strcmp(sprintf('%s%d',prefix,i),x.name),phys_names(phys_elecs));
    electrodes(i).nodes = unique(phys_names(phys_elecs(cur_elec)).nodes(:));
    electrodes(i).z_contact = z_contact;
end
