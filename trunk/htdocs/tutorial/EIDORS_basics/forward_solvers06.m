% Forward solvers $Id$

% Simulation Image
imgs= mk_image(mk_common_model('d2d1c',21));
imgs.fwd_model.electrode = imgs.fwd_model.electrode([15:21,1:8]);
nel = num_elecs(imgs);
imgs.fwd_solve.get_all_meas = 1;

% Output Image
imgo = rmfield(imgs,'elem_data');
imgo.calc_colours.ref_level = 0;

% Regular "current" stimulation between electrodes 6 and 10

stim.stim_pattern = zeros(nel,1);
stim.stim_pattern([6,10]) =  [10,-10];
stim.meas_pattern = speye(nel);
imgs.fwd_model.stimulation = stim;

vh = fwd_solve( imgs ); imgo.node_data = vh.volt(:,1);

subplot(221); show_fem(imgo,[0,1]); axis off;
