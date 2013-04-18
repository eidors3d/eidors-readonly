%Make an inverse model and extract forward model
imdl = mk_common_model('c2c0',16);
% 'b3cr',[16,3])
% imdl = mk_common_model('b3cr',[16,1]);
fmdl = imdl.fwd_model;
fmdl.nodes = fmdl.nodes*0.01;
% make point electrodes
% for i = 1:length(fmdl.electrode)
%     fmdl.electrode(i).nodes = fmdl.electrode(i).nodes(1);
% end


for i = 1:length(fmdl.stimulation);
   fmdl.stimulation(i).stim_pattern(fmdl.stimulation(i).stim_pattern > 0 ) = 0;
   fmdl.stimulation(i).stim_pattern(fmdl.stimulation(i).stim_pattern < 0 ) = 10;
   fmdl.stimulation(i).meas_pattern = sparse(diag(ones(16,1)));
   fmdl.stimulation(i).meas_pattern(i,i) = 0;
   fmdl.stimulation(i).meas_pattern(i,:) = [];
end
fmdl.meas_select = true(256,1);
fmdl.meas_select(1:17:256) = false;

%Default EIDORS solver
%Make image of unit conductivity
% img0 = mk_image(fmdl,0.15);
% img0.fwd_solve.get_all_meas = 1; %Internal voltage
% v0=fwd_solve(img0);
% v0e=v0.meas; v0all=v0.volt; 


%High-order EIDORS solver
%Change default eidors solvers
fmdl.solve = @fwd_solve_higher_order;
%fmdl.system_mat = @system_mat_higher_order; 
fmdl.system_mat = @calc_system_mat_opt;
% fmdl.jacobian = @jacobian_perturb;
fmdl.jacobian = @jacobian_adjoint_opt;
%fmdl.jacobian = @jacobian_adjoint_higher_order;

%Add element type
switch mdl_dim(fmdl)
   case 2
      fmdl.approx_type    = 'tri3';
   case 3
      fmdl.approx_type    = 'tet4'; % linear
end
%Make an image and get voltages using high order solver
fmdl.gnd_node = [];
img1 = mk_image(fmdl,15);
img1.fwd_solve.get_all_meas = 1; %Internal voltage
v1 = fwd_solve(img1); 
v1e=v1.meas; v1all=v1.volt;

imdl.fwd_model = fmdl;
%Plot electrode voltages and difference
%clf;subplot(211); plot([v0e,v1e]);
%legend('0','1'); xlim([1,256]);

v0 = fwd_solve(img1);
img1.elem_data(326) = 5;
imgl.gnd_node = [];
v1 = fwd_solve(img1);
rimg = inv_solve(imdl,v1,v0);
clf
% show_fem(rimg,[0 0 1]);
subplot(121)
show_fem(img1); view(2)
subplot(122)
switch mdl_dim(fmdl)
   case 3
      show_3d_slices(rimg); view(2)
   case 2
      show_fem(rimg);
end