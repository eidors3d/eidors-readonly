% Script to start EIDORS3D
% Set path and variables correctly

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: startup.m,v 1.13 2005-12-01 12:07:21 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms/aa_1996']);
addpath([HOMEDIR, '/algorithms/aa_2005']);
addpath([HOMEDIR, '/algorithms/ab_2002']);
%addpath([HOMEDIR, '/algorithms/at_2002']);
addpath([HOMEDIR, '/algorithms/ms_2005']);
%addpath([HOMEDIR, '/algorithms/mv_2001']);
addpath([HOMEDIR, '/algorithms/np_2003']);
addpath([HOMEDIR, '/interface']);
addpath([HOMEDIR, '/models/aa_1996']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
%addpath([HOMEDIR, '/tests']);

% test if eidors_var_id.cpp is a valid mexfile
if exist('eidors_var_id')~=3
  warning(sprintf([ ...
     'you do not have a compiled mex file eidors_var_id.\n' ...
     'Please compile it using: mex eidors_var_id.cpp\n'...
     'or if you have Matlab <6.5 under linux try:\n'...
     ' eval([''mex -v -f '' matlabroot ''/bin/cxxopts.sh eidors_var_id.cpp''])'
     ]));
end

clear HOMEDIR;

