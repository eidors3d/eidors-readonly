   ll=5; % length
   ww=1; % width
   hh=1; % height
   conduc= .13;  % conductivity in Ohm-meters
   current= 4;  % Amps
   z_contact= 9e-1;
   scale = .46;


   mdl= mk_grid_model([],0:ww,0:hh,0:ll);
   mdl.nodes= mdl.nodes*scale;

   mdl.gnd_node = 1;
   mdl.normalize_measurements = 0;
   
   elec_nodes= [1:(ww+1)*(hh+1)];
   elec(1).nodes= elec_nodes;
   elec(1).z_contact= z_contact;
   elec(2).nodes= num_nodes(mdl)-elec_nodes+1;
   elec(2).z_contact= z_contact;
   mdl.electrode= elec;
   show_fem(mdl)
   stim.stim_pattern= [-1;1]*current;
   stim.meas_pattern= [-1,1];
   mdl.stimulation= stim;


   
   mdl.solve = 'eidors_default';
   mdl.system_mat = 'eidors_default';

   img= mk_image(mdl, conduc);

   fsol= fwd_solve(img);
   fsol.meas
   
   %  show_fem(mdl);

   % analytical solution
   R = ll / ww / hh / scale/ conduc + 2*z_contact/scale^2;
   V= current*R