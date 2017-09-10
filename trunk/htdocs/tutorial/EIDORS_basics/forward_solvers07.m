% Forward solvers $Id$

% current stimulation between electrodes 6 and 10
stim.stim_pattern = zeros(nel,1);
stim.stim_pattern([6,10]) =  [10,-10];
% Voltage stimulation between electrodes 1,2 and 14,15
stim.volt_pattern = zeros(nel,1);
stim.volt_pattern([1,2,14,15]) =  [1,1,-1,-1]*5;
% Need to put NaN in the corresponding elements of stim_pattern
stim.stim_pattern([1,2,14,15]) =  NaN*ones(1,4);

imgs.fwd_model.stimulation = stim;
vh = fwd_solve( imgs ); imgo.node_data = vh.volt(:,1);

subplot(222); show_fem(imgo,[0,1]);
