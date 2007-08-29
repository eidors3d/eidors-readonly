% Script to start EIDORS
% Set path and variables correctly

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: startup.m,v 1.31 2007-08-29 09:16:48 aadler Exp $

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms']);
addpath([HOMEDIR, '/algorithms/a_adler']);
addpath([HOMEDIR, '/algorithms/a_borsic']);
addpath([HOMEDIR, '/algorithms/b_lionheart']);
%addpath([HOMEDIR, '/algorithms/m_vauhkonen']);
addpath([HOMEDIR, '/algorithms/n_polydorides']);
addpath([HOMEDIR, '/algorithms/d_stephenson']);
addpath([HOMEDIR, '/interface']);
addpath([HOMEDIR, '/models/a_adler']);
addpath([HOMEDIR, '/models/d_stephenson']);
addpath([HOMEDIR, '/models/s_murphy']);
addpath([HOMEDIR, '/meshing/netgen']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
%addpath([HOMEDIR, '/tests']);

% We need to add an architecture specific directory for mex files
archdir= '';
if exist('OCTAVE_VERSION')==5
   if findstr(computer,'-pc-');
      archdir= '/arch/octave/pc';
   end
else
   archdir= '/arch/matlab';
end
if ~isempty(archdir)
   addpath([HOMEDIR, archdir]);
end

% test if eidors_var_id.cpp is a valid mexfile
if exist('eidors_var_id')~=3
  warning(sprintf([ ...
     'you do not have a compiled mex file eidors_var_id.\n' ...
     'Please compile it using: mex eidors_var_id.cpp\n'...
     'or if you have Matlab <6.5 under linux try:\n'...
     ' eval([''mex -v -f '' matlabroot ''/bin/cxxopts.sh eidors_var_id.cpp''])'
     ]));
end

%prevent warnings in v7
if str2num(version('-release'))>=14
warning off MATLAB:symmmd:obsolete
end

% Setup defaults in calc_colours
calc_colours

% Set max cache size. Not completely sure about this
%  but 100MB should be available in most modern machines
eidors_cache('cache_size', 100e6 );


eidors_msg('Completed setting up of EIDORS Version %s', ...
    eidors_obj('eidors_version'),1);
eidors_msg('Parameter: cache_size=%d MB',eidors_cache('cache_size')/1e6,1);
eidors_msg('Parameter: mapped_colour=%d',calc_colours('mapped_colour'),1);
if calc_colours('greylev')>=0
   eidors_msg('Default background colour: white');
else
   eidors_msg('Default background colour: black');
end
eidors_msg('Architecture specific directory: %s',archdir,1);

clear HOMEDIR archdir;
