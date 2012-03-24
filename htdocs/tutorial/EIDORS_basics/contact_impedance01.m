Elec_width= 10; % 2 degrees - electrode width
params = [ 20,10,2]./[1000,1,100]; %d4 
ea = Elec_width/2 *(2*pi/360);
n_elec= 4;
for i=1:n_elec(1); 
  ai = (i-1)/n_elec(1) * 2*pi;
  elec_pts{i} = [sin(ai+ea),cos(ai+ea);sin(ai-ea),cos(ai-ea)];
end
fmdl= dm_2d_circ_pt_elecs( elec_pts, [], params);
subplot(221); show_fem(fmdl); axis image
print_convert contact_impedance01a.png
