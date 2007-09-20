% Build model
% $Id: build_single_plane01.m,v 1.1 2007-09-20 20:30:33 aadler Exp $

% Choose Parameters
electrodes_per_plane= 16;
number_of_planes= 1;
refine_electrodes= 10;
tank_radius= 15;
tank_height= 10;
electrode_width = 0.5;
electrode_height= 0.5;
rect_or_circ_electrode= 'C';


fname ='tank_mdl';

first_plane_starts= tank_height/(number_of_planes+1);
height_between_centres = first_plane_starts;

[fmdl,centres] = create_tank_mesh_ng( ...
   tank_radius, tank_height, ...
   rect_or_circ_electrode, ...
   log2(electrodes_per_plane), ...
   number_of_planes, ...
   first_plane_starts, ...
   height_between_centres, ...
   electrode_width, electrode_height, ...
   fname, refine_electrodes  );

% The msz file was created here can be reused later
msz_file= 'tank_mdl.msz';

% control mesh refinement: options are '-veryfine'; '-fine'; '';
finelevel= '';
if ~isempty(finelevel);
   call_netgen([fname,'.geo'],[fname,'.vol'],msz_file, finelevel);
end

% Create Homog
stim_pat= mk_stim_patterns(electrodes_per_plane, number_of_planes, ...
              '{ad}','{ad}',{'meas_current'});
[fmdl,mat_idxs]= ng_mk_fwd_model( ...
  [fname,'.vol'], centres, [], stim_pat);
n_elem= size(fmdl.elems,1);
elem_data= ones(n_elem,1);
img=eidors_obj('image','netgen_problem', ...
               'fwd_model',fmdl, 'elem_data', elem_data);
vh= fwd_solve( img);

show_fem( fmdl);
print -r100 -dpng build_single_plane01a.png;
