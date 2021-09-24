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
%  mat_indices:       cell array of element indices, per material

% Gmsh mesher for EIDORS was based on Netgen interface.
% (C) 2009 Bartosz Sawicki. License: GPL version 2 or version 3
% Modified by James Snyder, Bartlomiej Grychtol, Alistair Boyle

if ischar(vol_filename) && strcmp(vol_filename,'UNIT_TEST'); do_unit_test; return; end

if nargin<2;
   name = vol_filename;
end

if nargin<3 || isempty(eprefix);
   eprefix = 'electrode-';
end

if nargin < 4
    stim_pattern=[];
end

if nargin<5
   z_contact=0.01; % singular if z_contact=0
end


% Model Geometry
[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(vol_filename);
fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
                             stim_pattern, eprefix, z_contact, fc,phys_names);
[mat_indices,mat_names] = mk_mat_indices( mat_ind, phys_names);
if isempty(srf)
   fwd_mdl.boundary = find_boundary(fwd_mdl);
end

fwd_mdl.mat_idx = mat_indices;
if length(mat_names) > 0
    fwd_mdl.mat_names = mat_names;
end

% build fwd_model structure
function fwd_mdl= construct_fwd_model(srf,vtx,simp,bc, name, ...
    stim_pattern, eprefix, z_contact, fc, phys_names)
    mdl.nodes    = vtx;
    mdl.elems    = simp;
    mdl.boundary = srf;
    mdl.boundary_numbers=fc; % TODO this is not very useful without mapping to phys_names, like mat_idx/mat_names
    mdl.gnd_node = 1;
    mdl.name = name;

    % Model Stimulation
    if ~isempty(stim_pattern)
        mdl.stimulation= stim_pattern;
    end

    % Electrodes
    electrodes = find_elec(phys_names,eprefix,z_contact);
    if ~isempty(fieldnames(electrodes));
        mdl.electrode =     electrodes;
    end
    mdl.solve=          'eidors_default';
    mdl.jacobian=       'eidors_default';
    mdl.system_mat=     'eidors_default';

    fwd_mdl= eidors_obj('fwd_model', mdl);

% Output cell array of indices into each material type
%   array order is sorted by length of material type
function [mat_indices,mat_names]= mk_mat_indices(mat_ind,phys_names);
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
    mat_names = {};
    if length(phys_names) > 0
        phys_dims = max(cat(1,phys_names(:).dim));
        dim = max(phys_dims);
        mat_names = cell(1, l_idxs);
        for i = 1:l_idxs
            tag = sort_mi(find_mi(idxs(i)));
            idx = find(arrayfun(@(x) and((x.tag == tag), (x.dim == dim)), phys_names));
            assert(length(idx) == 1, 'missing physical name for tag');
            mat_names{i} = phys_names(idx).name;
        end
    end

% Assumes that electrodes are numbered starting at 1, with prefix provided
function electrodes = find_elec(phys_names,prefix,z_contact)
electrodes = struct();
phys_elecs = find(arrayfun(@(x)strncmp(x.name,prefix,length(prefix)),phys_names));
n_prefix = length(prefix);
for i = 1:length(phys_elecs)
    cur_elec = arrayfun(@(x) str2double(x.name((n_prefix+1):end)) == i,phys_names(phys_elecs));
    electrodes(i).nodes = unique(phys_names(phys_elecs(cur_elec)).nodes(:));
    electrodes(i).z_contact = z_contact;
end

function do_unit_test
DEBUG = 0; % enable if tests fail to help diagnose the problem

selfdir = fileparts(which('gmsh_read_mesh'));
vers = {'2.2', '4.0', '4.1'};
for ver = vers(:)'
    ver = ver{1};

    % Expected forward model:
    stim = 'asdf';
    elec.nodes = [];
    elec.z_contact = 0.11;
    fmdl2d.nodes = [0,0;1,0;1,1;0,1;0.5,0.5];
    fmdl2d.elems = [2,5,1;1,5,4;3,5,2;,4,5,3];
    fmdl2d.boundary = uint32([1,2;1,4;2,3;3,4]);
    fmdl2d.boundary_numbers = ones(4,1)*5;
    fmdl2d.gnd_node = 1;
    fmdl2d.name = 'test-name';
    fmdl2d.stimulation = stim;
    fmdl2d.solve = 'eidors_default';
    fmdl2d.jacobian = 'eidors_default';
    fmdl2d.system_mat = 'eidors_default';
    fmdl2d.electrode(1) = elec;
    fmdl2d.electrode(2) = elec;
    fmdl2d.electrode(1).nodes = [1,4]';
    fmdl2d.electrode(2).nodes = [2,3]';
    fmdl2d.type = 'fwd_model';
    fmdl2d.mat_idx = {{1:4}};
    fmdl2d.mat_names = {'main'};

    fmdl3d = fmdl2d;
    fmdl3d.electrode(1).nodes = [1,2,3,4,9]';
    fmdl3d.electrode(2).nodes = [5,6,7,8,10]';
    fmdl3d.mat_idx = {{1:24}};
    fmdl3d.boundary_numbers = [ones(4,1)*13; ones(4,1)*14];
    fmdl3d.boundary = [
         2     1     9
         1     3     9
         4     2     9
         3     4     9
         6    10     5
         5    10     7
         8    10     6
         7    10     8 ];
    fmdl3d.elems = [
        10    11    12    13
         9    12    14    11
        12    14    11    10
         9    12    11    13
         2     9     1    11
         1     9     3    14
        11    14     1     5
         4     9    12     3
         2     4     9    13
        12     3    14     7
         5    10    14     7
         7    10    12     8
        11    10     5     6
        12     4     8    13
        13     8    10     6
        13    11     2     6
        14     1     9    11
         9    12     3    14
        11     2     9    13
        12     9     4    13
         8    10    12    13
        14    11    10     5
        12    14    10     7
        13    10    11     6 ];
    fmdl3d.nodes = [
             0         0    1.0000
             0         0         0
             0    1.0000    1.0000
             0    1.0000         0
        1.0000         0    1.0000
        1.0000         0         0
        1.0000    1.0000    1.0000
        1.0000    1.0000         0
             0    0.5000    0.5000
        1.0000    0.5000    0.5000
        0.5000         0    0.5000
        0.5000    1.0000    0.5000
        0.5000    0.5000         0
        0.5000    0.5000    1.0000 ];

    % Test 2D and 3D parsing
    [fmdl,mat_ind] = gmsh_mk_fwd_model( ...
        fullfile(selfdir, ['box-' ver '.msh']), fmdl2d.name, 'elec#', stim, elec.z_contact );
    if DEBUG
        for x = fields(fmdl)'
            x = x{1};
            unit_test_cmp(['2d v' ver ' fmdl.' x],fmdl.(x),fmdl2d.(x))
        end
        fmdl
        fmdl2d
    end
    unit_test_cmp(['2d v' ver ' fmdl'],fmdl,fmdl2d)

    [fmdl,mat_ind] = gmsh_mk_fwd_model( ...
        fullfile(selfdir, ['cube-' ver '.msh']), fmdl3d.name, 'elec#', stim, elec.z_contact );
    if DEBUG
        for x = fields(fmdl)'
            x = x{1};
            unit_test_cmp(['3d v' ver ' fmdl.' x],fmdl.(x),fmdl3d.(x))
        end
        fmdl
        fmdl3d
    end
    unit_test_cmp(['3d v' ver ' fmdl'],fmdl,fmdl3d)
end
