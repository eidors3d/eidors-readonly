% Create 3D model of a tunnel $Id:$ 

N_elec = 16;
shape_str = ['solid incyl  = cylinder (0,0,0; 1,0,0; 1) -maxh=1.0; \n', ...
             'solid farcyl = cylinder (0,0,0; 1,0,0; 5) -maxh=5.0; \n' ...
             'solid pl1    =  plane(-5,0,0;-1,0,0);\n' ...
             'solid pl2    =  plane(5,0,0; 1,0,0);\n' ...
             'solid mainobj= pl1 and pl2 and farcyl and not incyl;\n'];
th= linspace(0,2*pi,N_elec+1)'; th(end)=[];
cth= cos(th); sth=sin(th); zth= zeros(size(th));
 elec_pos = [zth, cth, sth, zth cth, sth];
 elec_shape=[0.01];
 elec_obj = 'incyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
show_fem( fmdl );
view(90,0); print_convert tunnelsim01a.png
view(0,0); print_convert tunnelsim01b.png

crop_model([],  inline('x>=0.5','x','y','z'))
crop_model([],  inline('x<=-0.5','x','y','z'))
crop_model([],  inline('(y.^2+z.^2)>=1.3^2','x','y','z'))
view(-90,50); print_convert tunnelsim01c.png

