% Build model
% $Id: build_single_plane01.m,v 1.4 2007-09-21 19:07:47 aadler Exp $

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
msz_file= [fname, '.msz'];
stim_pat= mk_stim_patterns(electrodes_per_plane, number_of_planes, ...
              '{ad}','{ad}',{'meas_current'});

% control mesh refinement: options are '-veryfine'; '-fine'; '';
for finemodels= 0:2
  if     finemodels==0; finelevel= '';
  elseif finemodels==1; finelevel= '-fine';
  elseif finemodels==2; finelevel= '-veryfine';
  end

  if ~isempty(finelevel);
      call_netgen([fname,'.geo'],[fname,'.vol'],msz_file, finelevel);
  end

  [fmdl,mat_idxs]= ng_mk_fwd_model( [fname,'.vol'], centres, [], stim_pat);

  if     finemodels==0; ng_mdl_16x1_coarse= fmdl;
  elseif finemodels==1; ng_mdl_16x1_fine  = fmdl;
  elseif finemodels==2; ng_mdl_16x1_vfine = fmdl;
  end

  subplot(311);
  show_fem( fmdl); view(0,14);

  subplot(312);
  show_fem( fmdl); view(0,0);
  crop_model(gca, inline('y>0','x','y','z'))

  subplot(313);
  show_fem( fmdl); view(0,0);
  crop_model(gca, inline('y>-10','x','y','z'))
  set(gca,'Xlim',[-2,2],'Zlim',[-1,1]+first_plane_starts);

  print('-r100','-dpng', ...
        sprintf('build_single_plane01%c.png',96+finemodels));
end
   

save ng_mdl_16x1_coarse.mat  ng_mdl_16x1_coarse
save ng_mdl_16x1_fine.mat    ng_mdl_16x1_fine
save ng_mdl_16x1_vfine.mat   ng_mdl_16x1_vfine
