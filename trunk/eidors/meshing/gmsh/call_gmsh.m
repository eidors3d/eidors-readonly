function status= call_gmsh(geo_file, dim,extra)
% call_gmsh: call Gmsh to create a vol_file from a geo_file
% status= call_gmsh(geo_file, dim, extra)
%  staus = 0 -> success , negative -> failure
%
% geo_file = geometry file (input)
% dim      = model dimenstion (default: 2)
% extra    = additional string options for gmsh (see gmsh manual at
%            http://geuz.org/gmsh/doc/texinfo/gmsh.html#Command_002dline-options

% (C) 2009-2012 Bartosz Sawicki and Bartlomiej Grychtol
% License: GPL  version 2
% $Id$

% default to 2-D model
if nargin<2
    dim = 2;
end
if nargin<3
   extra = [];
end

cache_path = eidors_cache('cache_path');

if  exist('OCTAVE_VERSION') % FIXME
   islinux =1;
elseif ~strncmp(computer,'PC',2) % Don't know if we have isunix
   islinux =1;
else
   islinux =0;
end

% Gmsh executable filename
if  islinux
   gmsh_name = 'gmsh';
else
   gmsh_name = [cache_path,'/gmsh'];
end

while( 1 )
    ldpath='';
    if  ~exist('OCTAVE_VERSION') && ~isempty(strfind(system_dependent('getos'),'Linux'))
        islinux =1;
    else
        islinux =0;
    end
    
    status= system_cmd(sprintf( '%s "%s" -%d -v 2 %s',  gmsh_name, geo_file, dim, extra));
    
    if status==0; break; end
    
    if islinux
        disp(['It seems you are running Linux and Gmsh has not worked. ' ...
            'Check that it is installed and on the path. \n' ...
            'Perhaps LD_LIBRARY_PATH needs to be set?' ]);
        break;
    else
        fprintf([ ...
            'Gmsh call failed. Is Gmsh installed and on the search path?\n' ...
            'You are running under windows, I can attempt to create\n' ...
            'a batch file to access gmsh.\n' ...
            'Please enter the directory in which to find gmsh.\n' ...
            'If you dont have a copy, download it from' ...
            'http://www.geuz.org/gmsh/\n\n' ]);
        
        gmsh_path = input('gmsh_path? [or i=ignore, e=error] ','s');
        if strcmp(gmsh_name,'i'); break;end
        if strcmp(gmsh_name,'e'); error('user requested'),end;
        if exist( sprintf('%s/gmsh.exe',gmsh_path) , 'file' ) || ...
                exist( sprintf('%s/bin/gmsh.exe',gmsh_path) , 'file' )
            disp('Found gmsh');
            gmsh_exe = gmsh_path;
            if exist( sprintf('%s/bin/gmsh.exe',gmsh_path) , 'file' )
                gmsh_exe = [gmsh_path '/bin'];
            end
            
            
            fid= fopen([cache_path, '/gmsh.bat'],'w');
            fprintf(fid,'"%s/gmsh.exe" %%*\n', gmsh_exe);
            fclose(fid);
        end
    end
end
