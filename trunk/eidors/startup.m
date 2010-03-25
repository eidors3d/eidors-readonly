function startup
% Script to start EIDORS
% Set path and variables correctly

% NOTE: this is a function, so that we don't put variables into the
% workspace

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% CHECK WE HAVE THE RIGHT VERSION
% this is hard because matlab keeps on changing the format of the
% output of the version command. It used to be R13, now it is more
% like 2009a. Also the number of minor versions changes.

ver= version; ver(ver=='.')=' ';
ver = sscanf(ver,'%f'); ver=ver(:);
isoctave = exist('OCTAVE_VERSION')==5;

if ~isoctave
   if [1,1e-2]*ver(1:2) < 6.05
      warning(['EIDORS REQUIRES AT LEAST MATLAB V6.5.\n' ...
               'Several functions may not work with your version']);
   end
else
   if [1,1e-2,1e-4]*ver(1:3) < 3.0003
      warning(['EIDORS REQUIRES AT LEAST OCTAVE V3.0.3\n' ...
               'Several functions may not work with your version']);
   end
end

HOMEDIR=pwd;

addpath( HOMEDIR );
addpath([HOMEDIR, '/algorithms']);
addpath([HOMEDIR, '/algorithms/a_adler']);
addpath([HOMEDIR, '/algorithms/a_borsic']);
addpath([HOMEDIR, '/algorithms/b_lionheart']);
addpath([HOMEDIR, '/algorithms/c_gomez']);
addpath([HOMEDIR, '/algorithms/m_vauhkonen']);
addpath([HOMEDIR, '/algorithms/n_polydorides']);
addpath([HOMEDIR, '/algorithms/d_stephenson']);
addpath([HOMEDIR, '/algorithms/b_sawicki']);
addpath([HOMEDIR, '/interface']);
addpath([HOMEDIR, '/models/a_adler']);
addpath([HOMEDIR, '/models/d_stephenson']);
addpath([HOMEDIR, '/models/s_murphy']);
addpath([HOMEDIR, '/models/m_vauhkonen']);
addpath([HOMEDIR, '/meshing/netgen']);
addpath([HOMEDIR, '/meshing/distmesh']);
addpath([HOMEDIR, '/meshing/gmsh']);
addpath([HOMEDIR, '/sample_data']);
addpath([HOMEDIR, '/examples']);
addpath([HOMEDIR, '/tools']);
addpath([HOMEDIR, '/graphics_matlab']);
addpath([HOMEDIR, '/graphics_vtk']);
%addpath([HOMEDIR, '/tests']);

% We need to add an architecture specific directory for mex files
if isoctave 
   if findstr(computer,'x86_64-pc-');
      archdir= strcat('/arch/octave/',computer);
   elseif findstr(computer,'-pc-');
      archdir= '/arch/octave/pc';
   end
else
    % I don't know when matlab stopped using DLL as the extension
    % for WIN32 mex files. I'm guessing it's around 7.5
   if any(findstr(computer,'PCWIN')) && ( [1,1e-2]*ver(1:2) < 7.05 )
      archdir= '/arch/matlab/dll';
   else
      archdir= '/arch/matlab';
   end
end
addpath([HOMEDIR, archdir]);

% test if eidors_var_id.cpp is a valid mexfile
if exist('eidors_var_id')~=3
  if isoctave
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

% Setup defaults in calc_colours
calc_colours('greylev',-.001);          % background colour = white
calc_colours('sat_adj',.9);             % saturation of red and blue
calc_colours('window_range', .7);       % windowing of colours
calc_colours('backgnd',[.5,.5,.15]);    % background colour
calc_colours('mapped_colour',127);      % use 127*2+1 colourmap entries
calc_colours('ref_level','auto');       % auto set background colour
calc_colours('npoints',64);             % 64 raster points
calc_colours('clim',[]);                % no colour cropping
calc_colours('cb_shrink_move',[1,1,0]); % Don't shrink or move colorbar

% Set max cache size. Not completely sure about this
%  but 100MB should be available in most modern machines
eidors_cache('cache_size', 100e6 );
eidors_cache('boost_priority', 0 ); % set default priority

eidors_msg('Complete EIDORS (Ver: %s)', eidors_obj('eidors_version'),1);
eidors_msg('Parameter: cache_size=%d MB',eidors_cache('cache_size')/1e6,1);
eidors_msg('Parameter: mapped_colour=%d',calc_colours('mapped_colour'),1);
if calc_colours('greylev')>=0
   eidors_msg('Default background colour: black');
else
   eidors_msg('Default background colour: white');
end
eidors_msg('EIDORS mex folder: %s%s',HOMEDIR,archdir,1);

% check that the compiled mex file is newer than the source file
srcf = strcat(HOMEDIR,'/arch/eidors_var_id.cpp');
mexf = strcat(HOMEDIR,strcat(archdir,'/eidors_var_id.mex'));
if exist(srcf) == 2 & exist(mexf) == 2
  srcd=dir(srcf);
  mexd=dir(mexf);
  if datenum(srcd.date) > datenum(mexd.date)
     warning('eidors_var_id.mex file is older than source file and should be recompiled.');
  end
  clear srcd mexd;
end
clear srcf mexf;

% helpful messages
if ~exist('OCTAVE_VERSION');
   eidors_msg('New to EIDORS? Have a look at the <a href="http://eidors3d.sourceforge.net/tutorial/tutorial.shtml">Tutorials</a>.',1);
else 
   eidors_msg('New to EIDORS? Have a look at the Tutorials at http://eidors3d.sourceforge.net/tutorial/tutorial.shtml',1);
end
clear HOMEDIR archdir ver ans;

