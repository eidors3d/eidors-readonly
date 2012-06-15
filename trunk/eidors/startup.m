function startup( path_array )
% Script to start EIDORS
% Set path and variables correctly
% USAGE:
%   startup - setup basic eidors usage functions
%   startup( { dev directory paths })

% NOTE: this is a function, so that we don't put variables into the
% workspace

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% CHECK WE HAVE THE RIGHT VERSION
% this is hard because matlab keeps on changing the format of the
% output of the version command. It used to be R13, now it is more
% like 2009a. Also the number of minor versions changes.

ver= eidors_obj('interpreter_version');

if ver.isoctave
   if ver.ver < 3.002
      warning(['EIDORS REQUIRES AT LEAST OCTAVE V3.2.0\n' ...
               'Several functions may not work with your version']);
   end
else
   if ver.ver < 6.005
      warning(['EIDORS REQUIRES AT LEAST MATLAB V6.5.\n' ...
               'Several functions may not work with your version']);
   end
end

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/solvers']);
addpath([HOMEDIR, '/solvers/inverse']);
addpath([HOMEDIR, '/solvers/forward']);
addpath([HOMEDIR, '/algorithms']);
addpath([HOMEDIR, '/interface']);
addpath([HOMEDIR, '/models']);
addpath([HOMEDIR, '/meshing/netgen']);
addpath([HOMEDIR, '/meshing/distmesh']);
addpath([HOMEDIR, '/meshing/gmsh']);
addpath([HOMEDIR, '/meshing/stl']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/tools']);
addpath([HOMEDIR, '/graphics/matlab']);
addpath([HOMEDIR, '/graphics/vtk']);
addpath([HOMEDIR, '/external']);
addpath([HOMEDIR, '/deprecated']);
%addpath([HOMEDIR, '/tests']);

DEVDIR = [HOMEDIR(1:find(HOMEDIR == '/',1,'last')) '/dev'];
for i = 1:length(path_array)
    addpath([DEVDIR, '/', path_array{i}]);
end
% addpath([DEVDIR, '/a_adler']);
% addpath([DEVDIR, '/b_grychtol']);

% We need to add an architecture specific directory for mex files
if ver.isoctave 
   archdir= strcat('/arch/octave/',computer);
else
    % I don't know when matlab stopped using DLL as the extension
    % for WIN32 mex files. I'm guessing it's around 7.5
   if any(findstr(computer,'PCWIN')) && ( ver.ver < 7.005 )
      archdir= '/arch/matlab/dll';
   else
      archdir= '/arch/matlab';
   end
end
addpath([HOMEDIR, archdir]);

% test if eidors_var_id.cpp is a valid mexfile
if exist('eidors_var_id')~=3
  if ver.isoctave
    warning(sprintf([ ...
       'missing a required, pre-compiled mex file: eidors_var_id.\n' ...
       '  Please compile it using:\n'...
       '    cd ',HOMEDIR,'/arch\n'...
       '    mkoctfile -v --mex eidors_var_id.cpp\n'...
       '    mkdir -p ..',archdir,'\n'...
       '    mv *.mex ..',archdir
       ]));
   else
    warning(sprintf([ ...
       'missing a required, pre-compiled mex file: eidors_var_id.\n' ...
       '  Please compile it using:\n'...
       '     mex ',HOMEDIR,'/arch/eidors_var_id.cpp\n'...
       '  or if you have Matlab <6.5 under linux try:\n'...
       '    eval([''mex -v -f '' matlabroot ''/bin/cxxopts.sh eidors_var_id.cpp''])'
       ]));
   end
end

calc_colours('defaults'); % default calc_colours

% Set max cache size. Not completely sure about this
%  but 250MB should be available in most modern machines
eidors_cache('cache_size', 300e6 );
eidors_cache('boost_priority', 0 ); % set default priority

% Set default model cache location
mk_library_model('LIBRARY_PATH',[HOMEDIR, '/models/cache']);

eidors_msg('Installed EIDORS (Ver: %s)', eidors_obj('eidors_version'),1);
eidors_msg('Parameter: cache_size=%d MB',eidors_cache('cache_size')/1e6,1);
eidors_msg('Parameter: mapped_colour=%d',calc_colours('mapped_colour'),1);
if calc_colours('greylev')>=0
   eidors_msg('Default background colour: black');
else
   eidors_msg('Default background colour: white');
end
eidors_msg('EIDORS mex folder: %s%s',HOMEDIR,archdir,1);
eidors_msg('EIDORS model cache: %s', mk_library_model('LIBRARY_PATH'),1);

% check that the compiled mex file is newer than the source file
srcf = strcat(HOMEDIR,'/arch/eidors_var_id.cpp');
mexf = strcat(HOMEDIR,strcat(archdir,'/eidors_var_id.mex'));
if exist(srcf) == 2 && exist(mexf) == 2
  srcd=dir(srcf);
  mexd=dir(mexf);
  if srcd.datenum > mexd.datenum
     warning('eidors_var_id.mex file is older than source file and should be recompiled.');
  end
end

% helpful messages
% TODO: test if desktop is available
if ~ver.isoctave
   tutorials = '<a href="http://eidors3d.sf.net/tutorial/tutorial.shtml">Tutorials</a>';
else
   tutorials = 'Tutorials';
end
eidors_msg(['New to EIDORS? Have a look at the ',tutorials,'.'],1);

