% Script to start EIDORS3D
% Set path and variables correctly
% $Id: startup.m,v 1.7 2005-06-04 18:59:44 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms/np_2003']);
addpath([HOMEDIR, '/algorithms/aa_1996']);
addpath([HOMEDIR, '/models/aa_1996']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
addpath([HOMEDIR, '/tests']);

clear HOMEDIR;

