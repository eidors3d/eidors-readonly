function status= call_gmsh(geo_file)
% CALL_NETGEN: call netgen to create a vol_file from a geo_file
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
% $Id: call_netgen.m 1829 2009-06-08 21:30:00Z sawickib $
% (C) 2009 Bartosz Sawicki. Licensed under GPL V2


% Gmsh executable filename
gmsh_name = 'gmsh';

while( 1 )
   ldpath='';
   if  strfind(system_dependent('getos'),'Linux')
     islinux =1;
     s=version; ff= find(s=='.');
      if str2double(s(1:ff(2)-1))>=7
        %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks netgen    
          ldpath ='LD_LIBRARY_PATH=;';
      end      
   else
     islinux =0;
   end    

   status= system(sprintf( '%s %s %s -2 -v 2', ldpath, gmsh_name, geo_file));

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
%      gmsh_path = input('gmsh_path? ','s');
%      if exist( sprintf('%s/gmsh.exe',gmsh_path) , 'file' ) 
%         disp('Found gmsh.');
%       fid= fopen('ng.bat','w');
%       fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
%       fprintf(fid,'set TIX_LIBRARY=%s/lib/tix8.1\n', netgen_path);
%       fprintf(fid,'%s/netgen.exe %%*\n', netgen_path);
%       fclose(fid);
%       elseif exist( sprintf('%s/ng431.exe',netgen_path) , 'file' ) 
%          disp('Found netgen version 4.3.1');
%          fid= fopen('ng.bat','w');
%          fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
%          fprintf(fid,'set TIX_LIBRARY=%s/lib/tcl8.2\n', netgen_path);
%          fprintf(fid,'%s/ng431.exe %%*\n', netgen_path);
%          fclose(fid);
%       else
%          warning(['cannot find a version of netgen that I know about\n' ...
%                   'Install netgen 4.4 or 4.3.1 or check the path\n']);
%       end
        break;
   end
end
