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
clear;                      % Clear variables
%close all;                  % Close all windows
MyEIDORSStartup;            % Starts EIDORS

%eidors_cache('clear');      % Clear Cache

Script.tStart = tic;    % Remember start time

%% Forward Model - values in mm
Tank.Height = 236;
Tank.Radius = (242)/2;
Tank.maxMesh= 0;

Electrodes.NumberOf = 16;
Electrodes.ZPositions = 125;
Electrodes.InnerRadius1 = 0;
Electrodes.InnerRadius2 = 5;
Electrodes.OuterRadius1 = 10;
Electrodes.OuterRadius2 = 20;
Electrodes.maxh = 5;

fmdl = CreateTankWithCompoundElectrodes(Tank, Electrodes);

FEMimg = mk_image(fmdl, 1); 

% Show FEM forward model
figure(1); set(gcf, 'Name', 'FEM Forward Model'); show_fem(FEMimg, [0, 1, 0]); title('Forward Model');

%% Generate Stim Pattern
Pattern.NumberOfElectrodeRings = 1;
Pattern.AdjacentSkip = 7;
%Pattern.AdjacentSkip = 1;

%Pattern.RedundantMeasurementType = 'no_redundant';                       % 'do_redundant';
Pattern.RedundantMeasurementType = 'do_redundant';                      % 'do_redundant';
%Pattern.MeasureOnCurrentCarryingElectrodesType = 'no_meas_current';     % 'meas_current';
Pattern.MeasureOnCurrentCarryingElectrodesType = 'meas_current';         % 'meas_current';
Pattern.CurrentAmplitude = 5e-3;

fprintf('\n\nPattern Settings:\n---------\n');
fprintf('                       AdjacentSkip: %i\n', Pattern.AdjacentSkip);
fprintf('                       %s, %s\n', Pattern.RedundantMeasurementType, Pattern.MeasureOnCurrentCarryingElectrodesType);
fprintf('\n\n');

[fmdl.stimulation2, fmdl.meas_sel2] = mk_stim_patterns(Electrodes.NumberOf, 2*Pattern.NumberOfElectrodeRings, [0 Pattern.AdjacentSkip], [0 Pattern.AdjacentSkip], {Pattern.RedundantMeasurementType Pattern.MeasureOnCurrentCarryingElectrodesType}, Pattern.CurrentAmplitude);
[fmdl.stimulation, fmdl.meas_sel] = mk_stim_patterns(Electrodes.NumberOf, Pattern.NumberOfElectrodeRings, [0 Pattern.AdjacentSkip], [0 Pattern.AdjacentSkip], {Pattern.RedundantMeasurementType Pattern.MeasureOnCurrentCarryingElectrodesType}, Pattern.CurrentAmplitude);

for i=1:length(fmdl.stimulation)
    fmdl.stimulation(i).stim_pattern = full(fmdl.stimulation(i).stim_pattern);
    fmdl.stimulation(i).stim_pattern = [fmdl.stimulation(i).stim_pattern; zeros(16, 1)];
    fmdl.stimulation(i).stim_pattern = sparse(fmdl.stimulation(i).stim_pattern);
    
    fmdl.stimulation(i).meas_pattern = full(fmdl.stimulation(i).meas_pattern);
    fmdl.stimulation(i).meas_pattern = [zeros(16,16), fmdl.stimulation(i).meas_pattern];
    fmdl.stimulation(i).meas_pattern = [zeros(16,32); fmdl.stimulation(i).meas_pattern];
    fmdl.stimulation(i).meas_pattern = sparse(fmdl.stimulation(i).meas_pattern);
end
fmdl.stimulation3 = fmdl.stimulation2(1:16);

% Check im stim_pattern and stimulation is equal
for i=1:16
    if(all(fmdl.stimulation3(i).stim_pattern == fmdl.stimulation(i).stim_pattern) && all(fmdl.stimulation3(i).stimulation == fmdl.stimulation(i).stimulation))
        fprintf('%i ok\n',i);
    else
        fprintf('%i not ok\n',i);
    end;
end;

for i=1:16
    full(fmdl.stimulation(i).meas_pattern)
    full(fmdl.stimulation3(i).meas_pattern)
end;

ShowPattern(fmdl.stimulation, fmdl.meas_sel);

% DISPLAY_MEAS: display measurements on mesh
figure(2); display_meas(fmdl, 'ya', 0, 1);

%% Goodbye Message
Script.tElapsed=toc(Script.tStart);
fprintf('\nProcessed in %2.2f Seconds\n', Script.tElapsed);