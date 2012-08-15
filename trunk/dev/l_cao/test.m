%Make an inverse model and extract forward model
imdl = mk_common_model('c2C0',16);
fmdl = imdl.fwd_model;

% make point electrodes
for i = 1:length(fmdl.electrode)
    fmdl.electrode(i).nodes = fmdl.electrode(i).nodes(1);
end


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
img0 = mk_image(fmdl,0.15);
img0.fwd_solve.get_all_meas = 1; %Internal voltage
v0=fwd_solve(img0);
v0e=v0.meas; v0all=v0.volt; 


%High-order EIDORS solver
%Change default eidors solvers
fmdl.solve = @fwd_solve_higher_order;
% fmdl.system_mat = @system_mat_higher_order; 
fmdl.system_mat = @calc_system_mat_opt;
fmdl.jacobian = @jacobian_adjoint_opt;
imdl.fwd_model = fmdl;

%Add element type
fmdl.approx_type    = 'tri3'; % linear
%Make an image and get voltages using high order solver
img1 = mk_image(fmdl,0.15);
img1.fwd_solve.get_all_meas = 1; %Internal voltage
v1 = fwd_solve(img1); 
v1e=v1.meas; v1all=v1.volt;


%Plot electrode voltages and difference
clf;subplot(211); plot([v0e,v1e]);
legend('0','1'); xlim([1,256]);

v0 = fwd_solve(img1);
img1.elem_data(361) = 0.5;
v1 = fwd_solve(img1);
rimg = inv_solve(imdl,v1,v0);
clf
show_fem(rimg,[0 0 1]);