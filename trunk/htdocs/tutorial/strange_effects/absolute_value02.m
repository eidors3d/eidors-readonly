% Create Simulation Image
extra={'ball1', 'ball2','ball3',...
      ['solid ob = orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;' ...
       'solid ball1 = cylinder( 0.5, 0.2,0; 0.5, 0.2,1;0.2) and ob;' ...
       'solid ball2 = cylinder(-0.5, 0.2,0;-0.5, 0.2,1;0.2) and ob;' ...
       'solid ball3 = cylinder( 0.0,-0.5,0; 0.0,-0.5,1;0.2) and ob;']};

fmdl= ng_mk_cyl_models(0,[16],[0.1,0,0.05],extra); 
fmdl.stimulation = stim;
fmdl.meas_select = msel;
img = mk_image( fmdl, 1);
vh = fwd_solve(img); vh= vh.meas;

img.elem_data( fmdl.mat_idx{2} ) = 1.1;
img.elem_data( fmdl.mat_idx{3} ) = 0.9;
img.elem_data( fmdl.mat_idx{4} ) = 1.1;
vi = fwd_solve(img); vi= vi.meas;

subplot(221); show_fem(img);
print_convert absolute_value02a.png
