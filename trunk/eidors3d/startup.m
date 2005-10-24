% Script to start EIDORS3D
% Set path and variables correctly
% $Id: startup.m,v 1.8 2005-10-24 15:34:53 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms/np_2003']);
addpath([HOMEDIR, '/algorithms/aa_1996']);
addpath([HOMEDIR, '/algorithms/aa_2005']);
addpath([HOMEDIR, '/models/aa_1996']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
addpath([HOMEDIR, '/tests']);

% test if eidors_var_id.cpp is a valid mexfile
if exist('eidors_var_id')~=3
  warning(sprintf([ ...
     'you do not have a compiled mex file eidors_var_id.\n' ...
     'Please compile it using: mex eidors_var_id.cpp']));
end

clear HOMEDIR;

