clear;
close all;
clc;

% EIDORS has to be started before!

eidors_msg('log_level',3);

eidors_cache('cache_size', 4*1024*1024*1024);
eidors_cache('clear');
eidors_cache

fprintf('\nCache Test\n');

%% Generate Forward Model and Stim patterns

% Circular tank phantom in mm
TankHeight = 236;
TankRadius = 242/2;
cyl_shape = [TankHeight, TankRadius];        % Cylinder Shape = [Height, Radius, [maxsz]]

% Electrodes and position in mm
nrings = 1;                                  % Number of Electrode rings
nelec = 16;                                  % Number of Electrodes
ring_vert_pos = 125;                         % Position of the Electrode plances
elec_pos = [nelec, ring_vert_pos];
elec_shape = 5;                              % Electrode shape = [Width, Height, [maxsz]] ==> here circular electrodes 10 mm diameter

% Generate Model with Netgen
[fmdl, mat_idx] = ng_mk_cyl_models(cyl_shape, elec_pos, elec_shape);

%% Generate Stim Pattern
Pattern.NumberOfElectrodeRings = 1;
Pattern.AdjacentSkip = 7;
Pattern.RedundantMeasurementType = 'do_redundant';
Pattern.MeasureOnCurrentCarryingElectrodesType = 'meas_current';
Pattern.CurrentAmplitude = 5e-3;

% Create stim-pattern for 16 electrodes
[fmdl.stimulation, fmdl.meas_sel] = mk_stim_patterns(nelec, Pattern.NumberOfElectrodeRings, [0 Pattern.AdjacentSkip], [0 Pattern.AdjacentSkip], {Pattern.RedundantMeasurementType Pattern.MeasureOnCurrentCarryingElectrodesType}, Pattern.CurrentAmplitude);

% Add saved pattern to fmdl
fmdl.normalize_measurements = 0;

% Calculate homogenious forward solution
FEMimg = mk_image(fmdl, 1); 

% Show FEM forward model
%figure(); set(gcf, 'Name', 'FEM Forward Model'); show_fem(FEMimg, [0, 1, 0]); title('Forward Model');

%%
for i=1:5
    fprintf('\nIteration %i\n', i);
    eidors_cache('show_objs');
    vsim = fwd_solve(FEMimg);
end;