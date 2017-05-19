% based on http://eidors3d.sourceforge.net/tutorial/EIDORS_basics/contact_impedance.shtml
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

fmdl.stimulation(1).stim_pattern = [0;1;0;-1];
fmdl.stimulation(1).meas_pattern = [0;1;0;-1]';
fmdl.solve =      @fwd_solve_1st_order;
fmdl.system_mat = @system_mat_1st_order;
fmdl.electrode(1).z_contact = 0.01;

img = mk_image(fmdl,1);
img.fwd_solve.get_all_meas = 1;

imgc= img;
imgc.fwd_model.mdl_slice_mapper.npx = 128;
imgc.fwd_model.mdl_slice_mapper.npy = 200;
imgc.fwd_model.mdl_slice_mapper.level = [inf,inf,0];

img.fwd_model.electrode(1).z_contact=0.05;
img.fwd_model.stimulation(1).stim_pattern = [1;0;-1;0];
vh = fwd_solve(img);
imgc.fwd_model.mdl_slice_mapper.xpts = linspace(-0.25,0.25,200);
imgc.fwd_model.mdl_slice_mapper.ypts = linspace(0.8,1,100);
q = show_current(imgc,vh.volt);
hh=show_fem(imgc);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sx= linspace(-1,1,30); sy= sx*0;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;
axis([-.15,.15,0.85,1.02]);
axis off;
print('-dpdf','stream1.pdf');

img.fwd_model.electrode(1).z_contact=0.05;
img.fwd_model.stimulation(1).stim_pattern = [0;1;0;-1];
vh = fwd_solve(img);
imgc.fwd_model.mdl_slice_mapper.xpts = linspace(-0.25,0.25,200);
imgc.fwd_model.mdl_slice_mapper.ypts = linspace(0.8,1,100);
q = show_current(imgc,vh.volt);
hh=show_fem(imgc);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sy = linspace(.98,.8 ,20); sx= 0*sy - 0.15;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,-sx,sy); set(hh,'Linewidth',2);
hold off;
axis([-.15,.15,0.85,1.02]);
axis off;
print('-dpdf','stream2.pdf');
