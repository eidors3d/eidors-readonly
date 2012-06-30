function status= call_gmsh(geo_file, dim)
% call_gmsh: call Gmsh to create a vol_file from a geo_file
% status= call_netgen(geo_file, vol_file, msz_file, finelevel)
%  staus = 0 -> success , negative -> failure
%
% geo_file = geometry file (input)
% vol_file = FEM mesh file (output)
% msz_file = Meshsize file in netgen format
%
% Finelevel controls the fineness of the mesh
%   default is '' -> coarse
%   valid values are 'fine' or 'veryfine'
%
% $Id$
% (C) 2009 Bartosz Sawicki. Licensed under GPL V2

% default to 2-D model
if nargin<2
    dim = 2;
end

% Gmsh executable filename
gmsh_name = 'gmsh';

while( 1 )
   ldpath='';
   if  strfind(system_dependent('getos'),'Linux')
     islinux =1;
   else
     islinux =0;
   end    

   status= system_cmd(sprintf( '%s %s -%d -v 2',  gmsh_name, geo_file, dim));

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
       'http://www.geuz.org/gmsh/\n\n' ...
       'Note that you *MUST* use names without spaces. Thus\n' ...
       'instead of C:/Program Files/ write C:/Progra~1/\n\n' ]);
        break;
   end
end
