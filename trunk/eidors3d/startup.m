% Script to start EIDORS3D
% Set path and variables correctly
% $Id: startup.m,v 1.2 2004-07-02 14:07:12 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);

stderr=2; % file number of stderr
