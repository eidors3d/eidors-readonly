% Script to start EIDORS
% Set path and variables correctly

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: startup.m,v 1.18 2005-12-08 10:17:06 billlion Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms']);
addpath([HOMEDIR, '/algorithms/aa_1996']);
addpath([HOMEDIR, '/algorithms/aa_2005']);
addpath([HOMEDIR, '/algorithms/ab_2002']);
%addpath([HOMEDIR, '/algorithms/at_2002']);
addpath([HOMEDIR, '/algorithms/ms_2005']);
%addpath([HOMEDIR, '/algorithms/mv_2001']);
addpath([HOMEDIR, '/algorithms/np_2003']);
addpath([HOMEDIR, '/algorithms/ds_2005']);
addpath([HOMEDIR, '/interface']);
addpath([HOMEDIR, '/models/aa_1996']);
addpath([HOMEDIR, '/models/ds_2005']);
addpath([HOMEDIR, '/models/sm_2005']);
addpath([HOMEDIR, '/meshing/netgen']);
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

%prevent warnings in v7
if version('-release')>=14
warning off MATLAB:symmmd:obsolete
end
