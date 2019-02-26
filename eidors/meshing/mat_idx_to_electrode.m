function fmdl = mat_idx_to_electrode(fmdl, mat_idxes)
% MAT_IDX_TO_ELECTRODE: create electrodes from mat_idx values
% fmdl = mat_idx_to_electrode(fmdl, mat_idxes)
% fmdl: input and output fmdl
% Options: fmdl.mat_idx_to_electrode.z_contact  (z_contact value to use ... or default)
% mat_idxes = cell array of mat_idxes => convert electrode
%    example mat_idxes = {1:2, 5, [12,14]}

% (C) 2019 Andy Adler. License: GPL v2 or v3.
% $Id$

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

for i = 1:length(mat_idxes)
   fmdl = create_electrode_from_mat_idx(fmdl, mat_idxes{i});
end

fmdl = remove_unused_nodes(fmdl);

% adds electrode to the end
function fmdl = create_electrode_from_mat_idx(fmdl,nmat_idx);
   zc = 1e-5;
   try zc = fmdl.mat_idx_to_electrode.z_contact;
   end
   femobj = vertcat(fmdl.mat_idx{nmat_idx});
   for i=1:length(fmdl.mat_idx) 
      els = false(num_elems(fmdl),1);
      els(fmdl.mat_idx{i}) = 1;
      els(femobj) = [];
      fmdl.mat_idx{i} = find(els);
   end
   femobjnodes = fmdl.elems(femobj,:);
   fmdl.elems(femobj,:) = [];

   fmdl.boundary = find_boundary(fmdl);
   vt = intersect(femobjnodes,fmdl.boundary);
   elstr =  struct('nodes',vt(:)','z_contact',zc);
   if isfield(fmdl,'electrode')
     fmdl.electrode(end+1) = elstr;
   else
     fmdl.electrode(    1) = elstr;
   end


function do_unit_test
   extra={'ball_inside','ball_surface', [ ...
          'solid ball_inside  = sphere(-0.4,0,0.5;0.05);' ...
          'solid ball_surface = sphere( 0.4,0,1.0;0.05);' ...
          ]};
   fmdl= ng_mk_cyl_models(1,[8,.5],[.1],extra); 
   fmdl = mat_idx_to_electrode(fmdl, {2,3});
   unit_test_cmp('cyl_models',num_elecs(fmdl),10);

   show_fem(fmdl); view(3,12);
