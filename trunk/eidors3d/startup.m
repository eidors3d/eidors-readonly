% Script to start EIDORS3D
% Set path and variables correctly
% $Id: startup.m,v 1.6 2004-07-22 03:32:36 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms/np_2003']);
addpath([HOMEDIR, '/models/aa_1996']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
addpath([HOMEDIR, '/tests']);

clear HOMEDIR;

