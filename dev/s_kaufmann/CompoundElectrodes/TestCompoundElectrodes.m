% Module: TestCompoundElectrodes.m
%
% Test script for compound electrodes
%
% Known issues / Todo: -
%
% Input:    -
%
% Return:   -
%
% Written by Steffen Kaufmann <sk@steffen-kaufmann.com>, mainly based on 
% http://eidors3d.sourceforge.net/tutorial/netgen/netgen_gen_models.shtml
% Demo #17 from Andy Adler
% March 2014

%% Clear Workspace and Command Window, start time measurment & EIDORS
%clear all;                      % Clear variables
clc;
%close all;                  % Close all windows
%MyEIDORSStartup;            % Starts EIDORS

%eidors_cache('clear');      % Clear Cache

Script.tStart = tic;    % Remember start time

fprintf('CEV111 Image Reconstruction with EIDORS\n\n');

%% Parameters
Display.DisplayMeas = 0;
Display.ShowFEM = 1;
Display.ShowLegend = 1;
Display.SaveFigures = 1;

DoReconstruction = 1;
ChooseNF = 0;

UsePhase = 0;

bkgnd_conductivity =  0.65*1e-3;   % Tank bkgnd_conductivity set in order to align hom simulation to hom measurements
%bkgnd_conductivity =  1.155*1e-3;   % Tank bkgnd_conductivity set in order to align hom simulation to hom measurements

MeasurementFolder = 'C:\Repos\eidors\dev\s_kaufmann\CompoundElectrodes\Measurements\';
HomogeneousMeasFile = [MeasurementFolder 'SingleObject_Reference' '.mat'];
InhomogeneousMeasFile = [MeasurementFolder 'SingleObject_E1E2' '.mat'];

ExportDir = 'C:\Temp\';
ExportFile = [ExportDir 'CompoundElectrodes'];

prior = 'tikhonov';
%prior = 'noser';

%% Forward Model - values in mm

% Tank settings
Tank.Height = 236;
Tank.Radius = (242)/2;
Tank.maxMesh= 0;

% Electrode settings
Electrodes.NumberOf = 16;
Electrodes.ZPositions = 125;
Electrodes.InnerRadius1 = 0;
Electrodes.InnerRadius2 = 5;
Electrodes.OuterRadius1 = 10;
Electrodes.OuterRadius2 = 20;
Electrodes.maxh = 5;

% Scale Tank height in order to reduce the number of elements and to speed
% the solution up, with the trade off of beeing not accuratea anymore.
Scale = 1/5;
Tank.Height = Tank.Height * Scale;
Electrodes.ZPositions = Electrodes.ZPositions * Scale;

% The electrodes are numbered counter-clock wise
% Electrode 1..16 are the inner current electrodes
% Electrode 17..32 are the outer voltage electrodes
fmdl = CreateTankWithCompoundElectrodes(Tank, Electrodes);

% Reorder electrodes to be clock-wise arranged
fmdl.electrode([16:-1:1, 32:-1:17])= fmdl.electrode;

%% Generate Stim Pattern
Pattern.NumberOfElectrodeRings = 1;
Pattern.AdjacentSkip = 7;
%Pattern.AdjacentSkip = 1;

% For compound electrodes a full set can be accquired, because voltage and
% current electrodes are separated.
% Todo: Should be extracted out of the measurment files
%Pattern.RedundantMeasurementType = 'no_redundant';
Pattern.RedundantMeasurementType = 'do_redundant';
%Pattern.MeasureOnCurrentCarryingElectrodesType = 'no_meas_current';
Pattern.MeasureOnCurrentCarryingElectrodesType = 'meas_current';
Pattern.CurrentAmplitude = 5e-3;

% Create stim-pattern for 16 electrodes
[fmdl.stimulation, fmdl.meas_sel] = mk_stim_patterns(Electrodes.NumberOf, Pattern.NumberOfElectrodeRings, [0 Pattern.AdjacentSkip], [0 Pattern.AdjacentSkip], {Pattern.RedundantMeasurementType Pattern.MeasureOnCurrentCarryingElectrodesType}, Pattern.CurrentAmplitude);

% And scale them up for 32
for i=1:length(fmdl.stimulation)
    fmdl.stimulation(i).stim_pattern = full(fmdl.stimulation(i).stim_pattern);
    fmdl.stimulation(i).stim_pattern = [fmdl.stimulation(i).stim_pattern; zeros(16, 1)];
    fmdl.stimulation(i).stim_pattern = sparse(fmdl.stimulation(i).stim_pattern);
    
    fmdl.stimulation(i).meas_pattern = full(fmdl.stimulation(i).meas_pattern);
    fmdl.stimulation(i).meas_pattern = [zeros(16,16), fmdl.stimulation(i).meas_pattern];
    fmdl.stimulation(i).meas_pattern = [zeros(16,32); fmdl.stimulation(i).meas_pattern];
    fmdl.stimulation(i).meas_pattern = sparse(fmdl.stimulation(i).meas_pattern);
end

fprintf('\n\nPattern Settings:\n---------\n');
fprintf('                       AdjacentSkip: %i\n', Pattern.AdjacentSkip);
fprintf('                       %s, %s\n', Pattern.RedundantMeasurementType, Pattern.MeasureOnCurrentCarryingElectrodesType);
fprintf('\n\n');

ShowPattern(fmdl.stimulation, fmdl.meas_sel);

% Create Image of fwd mdl
fmdlIMG = mk_image(fmdl, bkgnd_conductivity, 'conductivity');

%% Simulate reference voltages according to the fwd mdl
simulation_data = fwd_solve(fmdlIMG);
%v_sim = simulation_data.meas(simulation_data.meas~=0);  % Remove 0 elements caused by not measure voltages on electrode 1:16 (current)
v_sim = simulation_data.meas;

%% Load Measurement Data
Data.Homogeneous = load(HomogeneousMeasFile);
Data.Inhomogeneous = load(InhomogeneousMeasFile);

%% Prepare Measurement Data
% Patch sign according to simulation and scale measured transfer impedance
% with I0 to get a voltage like EIDORS wants it.
% Todo: Should be done with Z_Cal in case the system was calibrated
Data.Homogeneous.v_eidors = (Data.Homogeneous.System.DAC.I0 * Data.Homogeneous.Z) .* sign(v_sim(simulation_data.meas~=0));
Data.Inhomogeneous.v_eidors = (Data.Inhomogeneous.System.DAC.I0 * Data.Inhomogeneous.Z) .* sign(v_sim(simulation_data.meas~=0));

% Scale up to 512 measurements 
v_hom = zeros(512, 1);
v_inhom = zeros(512, 1);
for i=1:16
    v_hom((17:32) + (32*(i-1))) = Data.Homogeneous.v_eidors( (1:16) .* i);
    v_inhom((17:32) +(32*(i-1))) = Data.Inhomogeneous.v_eidors( (1:16) .* i);
end;

%% Calculate Reciprocity
fprintf('\nCalculate Reciprocity...');
[RelativeReciprocityAccuracy_hom, AbsoluteReciprocityError_hom] = CalculateReciprocityAccuracy(Data.Homogeneous.System.Pattern, abs(Data.Homogeneous.v_eidors));
[RelativeReciprocityAccuracy_inhom, AbsoluteReciprocityError_inhom] = CalculateReciprocityAccuracy(Data.Inhomogeneous.System.Pattern, abs(Data.Homogeneous.v_eidors));
fprintf('done.\n');

%% Reconstruction
if DoReconstruction
    fprintf('\nStart reconstruction after %2.2f s....\n', toc(Script.tStart)); 

    % Create inverse model
    inv3d = eidors_obj('inv_model', 'EIT inverse');
    inv3d.reconst_type = 'difference';
    inv3d.jacobian_bkgnd.value = bkgnd_conductivity;        % What is this?
    inv3d.fwd_model = fmdl;
    %inv3d.fwd_model.np_fwd_solve.perm_sym = '{y}';         % What is this?
    inv3d.solve = @inv_solve_diff_GN_one_step;              % Is there a better algorithm?
    inv_model.inv_solve.calc_solution_error = 0;

    switch prior
        case 'tikhonov'
            % Tikhonov prior
            inv3d.R_prior = @prior_tikhonov;
            if ChooseNF
                inv3d = select_imdl(inv3d, {'Choose NF=1'});
            else
                inv3d.hyperparameter.value = 3.87e0;
            end;
        case 'noser'
            % Noser prior
            inv3d.R_prior = @prior_noser;
            inv3d.prior_noser.exponent=0.5; % std. value according to function header
            if ChooseNF
                inv3d = select_imdl(inv3d, {'Choose NF=1'});
            else
                inv3d.hyperparameter.value = 6.451171e0;
            end;
        case default
            error('Unkown Prior');
    end;
    fprintf(' uses: %s with a HP = %e...\n', prior, inv3d.hyperparameter.value);

    % Solve
    IMG_Solutions(1) = inv_solve(inv3d, v_hom, v_inhom);
    fprintf('Finished reconstruction after %2.2f s.\n', toc(Script.tStart));
end;

%% Ploting

%%% Display %%%
fprintf('\nDisplay\n');

%% Show FEM forward model
if (Display.ShowFEM)
    figure(1); set(gcf, 'Name', 'FEM Forward Model'); show_fem(fmdlIMG, [0, 1, 0]); title('Forward Model');
end;

%% Display Measurement Pattern
if (Display.DisplayMeas)
    figure(2); display_meas(fmdl, 'ya', 0, 1);
end;

%% Display RAW data
if (1)
    figure(3);  set(gcf, 'Name', 'U-Shapes');
    plot(1:length(v_sim(v_sim ~= 0)), 1e3*v_sim(v_sim ~= 0), '-o', 1:length(Data.Homogeneous.v_eidors), 1e3*Data.Homogeneous.v_eidors, '-x', 1:length(Data.Inhomogeneous.v_eidors), 1e3*Data.Inhomogeneous.v_eidors, '-x'); xlabel('No. Measurements'); ylabel('U / mV'); legend('Simulated_H_o_m', 'Measured_H_o_m', 'Measured_I_n_h_o_m');
end;

%% Display Reciprocity Error
if (length(v_sim) > 104)
    figure(4); set(gcf, 'Name', 'Reciprocity Accuracy');
    subplot(2,1,1); plot(1:length(Data.Homogeneous.v_eidors), 1e2*RelativeReciprocityAccuracy_hom, 'x-', 1:length(Data.Inhomogeneous.v_eidors), 1e2*RelativeReciprocityAccuracy_inhom, '-o'); xlabel('Measurement No.'); ylabel('Relative reciprocity accuracy / %'); legend('Measured_H_o_m', 'Measured_I_n_h_o_m');
    subplot(2,1,2); plot(1:length(Data.Homogeneous.v_eidors), 1e3*AbsoluteReciprocityError_hom, 'x-', 1:length(Data.Inhomogeneous.v_eidors), 1e3*AbsoluteReciprocityError_inhom, '-o'); xlabel('Measurement No.'); ylabel('Absolute reciprocity accuracy / mV'); legend('Measured_H_o_m', 'Measured_I_n_h_o_m');
end;

%% Display Reconstruction
if DoReconstruction
    for i=1:length(IMG_Solutions)
        IMG = IMG_Solutions(i);

        IMG.calc_colours.npoints = 256;
        IMG.calc_colours.cb_shrink_move = [1, 1, -0.1500]; % Leitfähigkeitsscala Balken (Width, Height, Position)

        Title = regexprep(IMG.name, '_', '\\_');

        figure();

        posn= [ inf,inf,Tank.Height,1; ...
                inf,inf,Electrodes.ZPositions-3,1; ...
                inf,inf,Electrodes.ZPositions,1 ; ...
                inf,inf,Electrodes.ZPositions+3,1; ...
                inf,inf,0,1];

        subplot(1, 2, 2); 
        show_slices(IMG, posn);

        if (Display.ShowLegend)
            ax = axes('position',[0,0,1,1], 'visible', 'off');
            text(0.82, 0.85,'Top');
            text(0.82, 0.69,'Electrode layer + 3 mm');
            text(0.82, 0.51,'Electrode layer');
            text(0.82, 0.35,'Electrode layer - 3 mm');
            text(0.82, 0.20,'Bottom');
        end;

        if UsePhase
            IMG.calc_colours.component = 'real';
            subplot(221); show_fem(IMG);title 'real conductivity change'

            IMG.calc_colours.component = 'imag'; 
            subplot(222); show_fem(IMG); title 'imag conductivity change'
        else
            subplot(1, 2, 1); show_fem(IMG, [1, 1]); title(Title);
        end;
    end;
    fprintf('done.\n');
end;

%% Save Figures
if (Display.SaveFigures)
    fprintf('Saving Figures....');
    h = get(0,'children');
    for i=1:length(h)
      saveas(h(i), [ExportFile num2str(i)], 'fig');
      saveas(h(i), [ExportFile num2str(i)], 'pdf');
    end;
    fprintf('done.\n');
end;

%% Goodbye Message
Script.tElapsed=toc(Script.tStart);
fprintf('\nProcessed in %2.2f Seconds\n', Script.tElapsed);