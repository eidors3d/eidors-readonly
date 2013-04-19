function eidors_startup( path_array )
% Script to start EIDORS
% Set path and variables correctly
% USAGE:
%   startup - setup basic eidors usage functions
%   startup( { dev directory paths })

% NOTE: this is a function, so that we don't put variables into the
% workspace

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin == 0
    path_array = {};
end

HOMEDIR=strrep(pwd,'\','/');
ver = version_check;
archdir = set_paths(HOMEDIR,ver, path_array);
eidors_cache('init');
set_defaults(HOMEDIR);
print_welcome(HOMEDIR,archdir, ver);

function set_defaults(HOMEDIR)
    % default functions
    eidors_default('set','fwd_solve','fwd_solve_1st_order');
    eidors_default('set','calc_system_mat','system_mat_1st_order');
    eidors_default('set','calc_jacobian','jacobian_adjoint');
    eidors_default('set','inv_solve','inv_solve_diff_GN_one_step');
    eidors_default('set','calc_RtR_prior','prior_laplace');
    eidors_default('set','calc_R_prior','prior_laplace');
    eidors_default('set','fwd_solve','fwd_solve_1st_order');
    %models are NOT normalized by default
    eidors_default('set','mdl_normalize',@(x) 0); 
    


    calc_colours('defaults'); % default calc_colours

    % Set max cache size. Not completely sure about this
    %  but 250MB should be available in most modern machines
    eidors_cache('cache_size', 300e6 );
    eidors_cache('boost_priority', 0 ); % set default priority

    % Set default model cache location
    mk_library_model('LIBRARY_PATH',[HOMEDIR, '/models/cache']);


function ver = version_check
    ver= eidors_obj('interpreter_version');

    if ver.isoctave
        if ver.ver < 3.002
            warning(['EIDORS REQUIRES AT LEAST OCTAVE V3.2.0\n' ...
                'Several functions may not work with your version']);
        end
    else
        if ver.ver < 7.000
            warning(['EIDORS REQUIRES AT LEAST MATLAB V7.0.\n' ...
                'Several functions may not work with your version']);
        end
    end

function archdir = set_paths(HOMEDIR, ver,path_array)

    addpath( HOMEDIR );
    addpath([HOMEDIR, '/solvers']);
    addpath([HOMEDIR, '/solvers/inverse']);
    addpath([HOMEDIR, '/solvers/forward']);
    addpath([HOMEDIR, '/algorithms']);
    addpath([HOMEDIR, '/interface']);
    addpath([HOMEDIR, '/models']);
    addpath([HOMEDIR, '/meshing']);
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
    addpath([HOMEDIR, '/overloads']);

    %addpath([HOMEDIR, '/tests']);

    DEVDIR = [HOMEDIR(1:find(HOMEDIR == filesep,1,'last')-1) '/dev'];
    for i = 1:length(path_array)
        p = genpath([DEVDIR, '/', path_array{i}]);
        addpath(p);
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
    fname = [HOMEDIR, archdir, '/eidors_var_id.', mexext];
    
    if ~exist(fname, 'file')
       warning('missing a required, pre-compiled mex file: eidors_var_id');
       compile_mex(HOMEDIR,archdir,ver);
    end

    % check that the compiled mex file is newer than the source file
    srcf = strcat(HOMEDIR,'/arch/eidors_var_id.cpp');
    mexf = strcat(fname);
    if exist(srcf) == 2 && exist(mexf) == 3
        srcd=dir(srcf);
        mexd=dir(mexf);
        if srcd.datenum > mexd.datenum
           if ver.isoctave
              warning(sprintf([ ...
                 'eidors_var_id.mex file is older than source file and should be recompiled.\n' ...
                 '  Please compile it using:\n'...
                 '    cd ',HOMEDIR,'/arch\n'...
                 '    mkoctfile -v --mex eidors_var_id.cpp\n'...
                 '    mkdir -p ..',archdir,'\n'...
                 '    mv *.mex ..',archdir
                 ]));
           else
              warning(sprintf([ ...
                 'eidors_var_id.mex file is older than source file and should be recompiled.\n' ...
                 '  Please compile it using:\n'...
                 '     cd ',HOMEDIR,'/arch\n'...
                 '     mex ',HOMEDIR,'/arch/eidors_var_id.cpp\n'...
                 '     mv *.mex* ',HOMEDIR,'/arch/matlab\n' ...
                 'If you have a 64 bit machine, please use "mex -largeArrayDims ..."\n' ...
                 ]));
           end
        end
    end

function compile_mex(HOMEDIR,archdir, ver)
    eidors_msg('Attempting to compile eidors_var_id',2);
    c = computer;
    flags = [];

    if ver.isoctave
         curdir = cd;
         cd(sprintf('%s/arch',HOMEDIR));
         mkoctfile -v --mex eidors_var_id.cpp
         eval(sprintf('mkdir -p ..%s',archdir));
         movefile(sprintf('%s/arch/*.mex',HOMEDIR), ...
                  sprintf('%s%s',HOMEDIR,archdir));
         cd(curdir)
         return
    end
    
    if strcmp(c(end-1:end),'64')
       flags = '-largeArrayDims';
    end  
    cmd = sprintf('mex %s %s/arch/eidors_var_id.cpp', flags, HOMEDIR);
    eval(cmd);
    movefile(sprintf('%s/*.mex*',HOMEDIR), ...
           sprintf('%s%s',HOMEDIR,archdir));

function print_welcome(HOMEDIR,archdir,ver)
    eidors_ver = eidors_obj('eidors_version');
    if eidors_ver(end) == '+' % post release version
       % THIS IS HORRIBLE, HORRIBLE CRAP IN SVN. LOTS OF USERS WANT GlobalRev
       % BUT THE ARROGANT SVN AUTHORS REFUSE TO PROVIDE IT!!!!
       [status, result] = system('svnversion');
       if status==0;
          eidors_ver = [eidors_ver, ' SVN_ID=', result(1:end-1)];
       end
    end
    eidors_msg('Installed EIDORS (Ver: %s)', eidors_ver,1);

    eidors_msg('Parameter: cache_size=%d MB',eidors_cache('cache_size')/1e6,1);
    eidors_msg('Parameter: mapped_colour=%d',calc_colours('mapped_colour'),1);
    if calc_colours('greylev')>=0
        eidors_msg('Default background colour: black',1);
    else
        eidors_msg('Default background colour: white',1);
    end
    eidors_msg('EIDORS mex folder: %s%s',HOMEDIR,archdir,1);
    eidors_msg('EIDORS model cache: %s', mk_library_model('LIBRARY_PATH'),1);


    % helpful messages
    % TODO: test if desktop is available
    if ver.isoctave
        canwritehtml=0;
    else try
            mf = com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame;
            if isempty(mf)
                canwritehtml=0;
            else
                canwritehtml=1;
            end
        catch
            canwritehtml=0;
        end
    end
    if canwritehtml
        tutorials = '<a href="http://eidors3d.sf.net/tutorial/tutorial.shtml">Tutorials</a>';
    else
        tutorials = 'Tutorials';
    end
    eidors_msg(['New to EIDORS? Have a look at the ',tutorials,'.'],1);

