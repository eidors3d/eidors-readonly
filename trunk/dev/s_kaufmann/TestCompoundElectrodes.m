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
clc;
%close all;                  % Close all windows
MyEIDORSStartup;            % Starts EIDORS

%eidors_cache('clear');      % Clear Cache

Script.tStart = tic;    % Remember start time

%% Parameters
Display.DisplayMeas = 0;
Display.ShowFEM = 1;

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

% The electrodes are numbered counter-clock wise
% Electrode 1..16 are the inner current electrodes
% Electrode 17..32 are the outer voltage electrodes
fmdl = CreateTankWithCompoundElectrodes(Tank, Electrodes);

% Reorder electrodes to be clock-wise arranged
fmdl.electrode([16:-1:1, 32:-1:17])= fmdl.electrode;

% Create image
FEMimg = mk_image(fmdl, 1); 

% Show FEM forward model
if (Display.ShowFEM)
    figure(1); set(gcf, 'Name', 'FEM Forward Model'); show_fem(FEMimg, [0, 1, 0]); title('Forward Model');
end;

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

% Create stim-pattern for 16 electrodes
[fmdl.stimulation, fmdl.meas_sel] = mk_stim_patterns(Electrodes.NumberOf, Pattern.NumberOfElectrodeRings, [0 Pattern.AdjacentSkip], [0 Pattern.AdjacentSkip], {Pattern.RedundantMeasurementType Pattern.MeasureOnCurrentCarryingElectrodesType}, Pattern.CurrentAmplitude);

% And scale them up for 30
for i=1:length(fmdl.stimulation)
    fmdl.stimulation(i).stim_pattern = full(fmdl.stimulation(i).stim_pattern);
    fmdl.stimulation(i).stim_pattern = [fmdl.stimulation(i).stim_pattern; zeros(16, 1)];
    fmdl.stimulation(i).stim_pattern = sparse(fmdl.stimulation(i).stim_pattern);
    
    fmdl.stimulation(i).meas_pattern = full(fmdl.stimulation(i).meas_pattern);
    fmdl.stimulation(i).meas_pattern = [zeros(16,16), fmdl.stimulation(i).meas_pattern];
    fmdl.stimulation(i).meas_pattern = [zeros(16,32); fmdl.stimulation(i).meas_pattern];
    fmdl.stimulation(i).meas_pattern = sparse(fmdl.stimulation(i).meas_pattern);
end

%% Display
% Plot Pattern
ShowPattern(fmdl.stimulation, fmdl.meas_sel);

% Display Measurement Pattern
if (Display.DisplayMeas)
    figure(2); display_meas(fmdl, 'ya', 0, 1);
end;

%% Goodbye Message
Script.tElapsed=toc(Script.tStart);
fprintf('\nProcessed in %2.2f Seconds\n', Script.tElapsed);