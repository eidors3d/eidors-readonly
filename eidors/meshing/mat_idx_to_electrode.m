function [fmdl,rm_elems] = mat_idx_to_electrode(fmdl, mat_idxes)
% MAT_IDX_TO_ELECTRODE: create electrodes from mat_idx values
% fmdl = mat_idx_to_electrode(fmdl, mat_idxes)
% fmdl: input and output fmdl
% Options: fmdl.mat_idx_to_electrode.z_contact  (z_contact value to use ... or default)
% mat_idxes = cell array of mat_idxes => convert electrode
%    example mat_idxes = {1:2, 5, [12,14]}
%
% by default, adds a 'faces'-type electrode (i.e. an
%    electrode defined in terms of its faces)
% Parameter
%   fmdl.mat_idx_to_electrode.nodes_electrode = true
%      Add a 'nodes'-type electrode
%
% To work with an image object, use:
%  [img.fwd_model,rm_elems] = mat_idx_to_electrode(img.fwd_model,{mat_idxes});
%  img.elem_data(rm_elems) = [];
%
% Notes:
%   mat_idx regions cannot be directly touching each other.

% (C) 2019 Andy Adler. License: GPL v2 or v3.
% $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

elec_faces = true;
try 
   elec_faces = ~fmdl.mat_idx_to_electrode.nodes_electrode;
end

rm_elems = [];
for i = 1:length(mat_idxes)
   if elec_faces 
      [fmdl,rm_elemi] = create_electrode_faces_from_mat_idx(fmdl, mat_idxes{i});
   else
      [fmdl,rm_elemi] = create_electrode_nodes_from_mat_idx(fmdl, mat_idxes{i});
   end
   rm_elems = union(rm_elems,rm_elemi); 
end

fmdl = remove_unused_nodes(fmdl);
fmdl = remove_unused_boundary(fmdl);
if any(fmdl.boundary(:) == 0)
   eidors_msg('WARNING: PROBLEM WITH BOUNDARY',1)
   keyboard
end

% This code adds the new electrode as a elec(i).nodes
% adds electrode to the end
% femobj is removed elems
function [fmdl,femobj] = create_electrode_nodes_from_mat_idx(fmdl,nmat_idx);
   zc = 1e-5;
   try zc = fmdl.mat_idx_to_electrode.z_contact;
   end
   femobj = vertcat(fmdl.mat_idx{nmat_idx});
   faces = calc_elec_faces(fmdl.elems,femobj);
   face1=faces;

   femobjnodes = fmdl.elems(femobj,:);
   [fmdl,faces] = rm_elems( fmdl, femobj, faces);

   % Add to the boundary with the boundary of the new electrode
   fmdl.boundary = unique([fmdl.boundary;faces],'rows');

   vt = find_bdynodes(fmdl,femobjnodes);
   elstr =  struct('nodes',vt(:)','z_contact',zc);
   fmdl = add_elec(fmdl, elstr);

function fmdl = remove_unused_boundary(fmdl);
   elems= fmdl.elems;
   fmdl.boundary = unique(sort(fmdl.boundary,2),'rows');
   switch mdl_dim(fmdl)
      case 2; selem = [sort(elems(:,[1,2]),2);
                       sort(elems(:,[1,3]),2);
                       sort(elems(:,[2,3]),2)];
      case 3; selem = [sort(elems(:,[1,2,3]),2);
                       sort(elems(:,[1,2,4]),2);
                       sort(elems(:,[1,3,4]),2);
                       sort(elems(:,[2,3,4]),2)];
      otherwise; error('Dims must be 2 or 3');
   end
   [~,idx] = setdiff(fmdl.boundary, selem,'rows');
%  disp( fmdl.boundary(idx,:) )
   fmdl.boundary(idx,:) = [];

function [fmdl,faces] = rm_elems( fmdl, femobj, faces);
   % fix the mat_idx object, since we remove femobj
   for i=1:length(fmdl.mat_idx) 
      els = false(num_elems(fmdl),1);
      els(fmdl.mat_idx{i}) = true;
      els(femobj) = [];
      fmdl.mat_idx{i} = find(els);
   end
   fmdl.elems(femobj,:) = [];

% Taken care of by remove_unused_boundary
%  % remove faces that are no longer on the boundary
%  if nargin>2
%     inels = reshape(... 
%          ismember(faces(:),fmdl.elems(:)), size(faces));
%     faces(any(~inels,2),:) = [];
%     sfaces = sort(faces,2);
%     selem1 = sort(fmdl.elems(:,[1,2,3]),2);
%     selem2 = sort(fmdl.elems(:,[1,2,4]),2);
%     selem3 = sort(fmdl.elems(:,[1,3,4]),2);
%     selem4 = sort(fmdl.elems(:,[2,3,4]),2);
%     [~,idx] = setdiff(sfaces,[selem1;selem2;selem3;selem4],'rows');
%     faces(idx,:) = [];
%  end

% This code adds the new electrode as a elec(i).faces
% adds electrode to the end
% femobj is removed elems
function [fmdl,femobj] = create_electrode_faces_from_mat_idx(fmdl,nmat_idx);
   zc = 1e-5;
   try zc = fmdl.mat_idx_to_electrode.z_contact;
   end
   femobj = vertcat(fmdl.mat_idx{nmat_idx});
   faces = calc_elec_faces(fmdl.elems,femobj);
   [fmdl,faces] = rm_elems( fmdl, femobj, faces);
   fmdl.boundary = unique([fmdl.boundary;faces],'rows');

   elstr =  struct('nodes',[],'z_contact',zc,'faces',faces);

   fmdl = add_elec(fmdl, elstr);

function fmdl = add_elec(fmdl, elstr)
   if isfield(fmdl,'electrode')
%    Stupid matlab forces you to add this
     if isfield(elstr,'faces') && ~isfield(fmdl.electrode,'faces');
        [fmdl.electrode(:).faces] = deal([]);
     end
     fmdl.electrode(end+1) = elstr;
   else
     fmdl.electrode(    1) = elstr;
   end

function faces = calc_elec_faces(elems,femobj);
   thisels = elems(femobj,:);
   switch size(thisels,2)  % dimentions
     case 2; % 1D
       allfaces = thisels;
     case 3; % 2D
       allfaces = [thisels(:,[1,2]);
                   thisels(:,[2,3]);
                   thisels(:,[1,3])];
     case 4; % 3D
       allfaces = [thisels(:,[1,2,3]);
                   thisels(:,[1,2,4]);
                   thisels(:,[1,3,4]);
                   thisels(:,[2,3,4])];
     otherwise;
       error 'Dimensions of elems';
   end
   allfaces = sort(allfaces,2);
   % remove all faces which exist more than once
   [faces,ia,ib] = unique(allfaces,'rows');
   exist_faces = sparse(1,ib,1);
   faces(exist_faces>1,:) = [];

function vt = find_bdynodes(fmdl,femobjnodes)
% Slow way
%  vt = intersect(femobjnodes,find_boundary(fmdl));
% Fast: pre-process 
   usenodes = reshape( ismember( ...
      fmdl.elems, femobjnodes), size(fmdl.elems));
   fmdl.elems(~any(usenodes,2),:) = [];
   vt = intersect(femobjnodes,find_boundary(fmdl));
      

function do_unit_test
   clf; subplot(221);
   do_unit_test_2d(true)
   do_unit_test_2d(false)
   subplot(223);
   do_unit_test_3d_netgen
   do_unit_test_3d_netgen2

function do_unit_test_2d(nodes_electrode)
   fmdl = getfield(mk_common_model('a2c2',1),'fwd_model');
   fmdl.mat_idx_to_electrode.nodes_electrode = nodes_electrode;
   fmdl.electrode(1).nodes = 26:41;
   fmdl.mat_idx{1} = [1:4];
   fmdl= mat_idx_to_electrode(fmdl, {1});
   unit_test_cmp('a2c2-01',num_elems(fmdl),64-4);
   if nodes_electrode;
      unit_test_cmp('a2c2-02',length(fmdl.electrode(2).nodes),4);
   else
      unit_test_cmp('a2c2-02a',length(fmdl.electrode(2).nodes),0);
      unit_test_cmp('a2c2-02b',size(fmdl.electrode(2).faces),[4,2]);
   end

   fmdl.stimulation = stim_meas_list([1,2,1,2]);
   img = mk_image(fmdl,1);
   img.fwd_solve.get_all_meas = true;
   vh = fwd_solve(img);

   img = rmfield(img,'elem_data');
   img.node_data = vh.volt;
   
   img.fwd_model = rmfield(img.fwd_model,'electrode');
   show_fem(img,[0,0,2]);
%  plot(sqrt(sum(fmdl.nodes.^2,2)),vh.volt,'*');
   unit_test_cmp('a2c2-03', std(vh.volt(13:24)),0,0.002);
   vd = mean(vh.volt(5:12)) - mean(vh.volt(13:24));
   unit_test_cmp('a2c2-03', vd,0.065251,1e-6);



function do_unit_test_3d_netgen
   extra={'ball_inside','ball_surface', [ ...
          'solid ball_inside  = sphere(-0.4,0,0.5;0.05);' ...
          'solid ball_surface = sphere( 0.4,0,1.0;0.05) -maxh=.005;' ...
          ]};
   fmdl= ng_mk_cyl_models(1,[8,.5],[.1],extra); 
   fmdl.mat_idx_to_electrode.nodes_electrode = true;

   fmd_= mat_idx_to_electrode(fmdl, {2,3});
   unit_test_cmp('cyl_models 01',num_elecs(fmd_),10);

   fmd_= mat_idx_to_electrode(fmdl, {2:3});
   unit_test_cmp('cyl_models 02',num_elecs(fmd_),9);

   fmd_= mat_idx_to_electrode(fmdl, {2});
   unit_test_cmp('cyl_models 03',num_elecs(fmd_),9);

   img = mk_image(fmd_,1);
   img.elem_data(fmd_.mat_idx{3})= 1.1;
   img.calc_colours.ref_level = 1;
   unit_test_cmp('cyl_models 04',num_elecs(img),9);

   img = mk_image(fmdl,1);
   img.elem_data(vertcat(fmdl.mat_idx{2:3}))= 1.1;
   [img.fwd_model,rm_elems]= mat_idx_to_electrode( ...
        img.fwd_model, {2});
   img.elem_data(rm_elems) = [];
   unit_test_cmp('cyl_models 05',num_elecs(img),9);
   vol = get_elem_volume(img);
   vvs = sum(vol(find(img.elem_data == 1.1)));
   unit_test_cmp('cyl_models 05a',vvs,4/3*pi*.05^3,1e-5);
  


   img.calc_colours.ref_level = 1;
   show_fem(img); view(3,12);


function do_unit_test_3d_netgen2
   shape_str = ['solid sqelec = orthobrick(-1,-1,-1;1,1,0); tlo sqelec;' ...
                'solid mainobj= orthobrick(-5,-5,-5;5,5,0) and not sqelec;'];
   fmdl = ng_mk_gen_models(shape_str, [],[],'');
   fmdl.mat_idx_to_electrode.nodes_electrode = true;
   fmdl = mat_idx_to_electrode(fmdl,{1});
   elimnodes = [2 35 36; 4 34 36; 6 31 35; 8 31 34; 31 34 35; 34 35 36];
   sd = intersect(elimnodes,sort(fmdl.boundary,2),'rows');
   unit_test_cmp('cut boundary', size(sd,1),0);

