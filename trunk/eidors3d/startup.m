% Script to start EIDORS3D
% Set path and variables correctly
% $Id: startup.m,v 1.3 2004-07-17 14:26:13 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
addpath([HOMEDIR, '/tests']);

clear HOMEDIR;
