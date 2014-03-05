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
close all;                  % Close all windows
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

%% Goodbye Message
Script.tElapsed=toc(Script.tStart);
fprintf('\nProcessed in %2.2f Seconds\n', Script.tElapsed);