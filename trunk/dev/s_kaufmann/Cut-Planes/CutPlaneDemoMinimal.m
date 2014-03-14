% Module: CutPlaneDemo.m
%
% Displays some stuff
%
% Known issues / Todo: Requires started EIDORS
%
% Input:    -
%
% Return:   -
%
% Written by Steffen Kaufmann <sk@steffen-kaufmann.com>
% March 2014

%% Clear Workspace
close all;
clc;
clear;

% Electrode Position
Electrodes.ZPositions = 130;

% Load reconstructed Image
load('C:\Repos\EIDORS\dev\s_kaufmann\Cut-Planes\CutPlanes.mat');

% Display 
figure(); set(gcf, 'Name', 'Show_3d_slices - 2');
hold on;

th = linspace(0,2*pi, 16+1)'; th(end)=[];
on_elecs = ones(16, 1);
el_th = []; 
el_z  = []; 
el_th = [el_th; th];
el_z  = [el_z ; on_elecs*Electrodes.ZPositions+5];

elec_pos = [242/2*cos(el_th), 242/2*sin(el_th), el_z];    

% Draw empty FEM with electrode numbers
show_fem(IMG_Solutions(1).fwd_model, [0, 1, 0]);
view(3);

for i=1:size(elec_pos, 1)
    slc = mdl_slice_mesher(IMG_Solutions(1), elec_pos(i,:));
    slc.calc_colours = IMG_Solutions(1).calc_colours;
    show_fem(slc);
    fprintf('%i. Press any key to continue...\n', i);  pause;
end;
hold off;
