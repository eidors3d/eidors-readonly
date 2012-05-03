function status= call_netgen(geo_file, vol_file, msz_file, finelevel)
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
% $Id$
% (C) 2006 Andy Adler. Licensed under GPL V2

if nargin<3
   msz_file= '';
end

if ~isempty(msz_file)
    eidors_msg('call_netgen: Warning. Using an *.msz file. This often fails.');
end

if nargin<4
%  finelevel= '-veryfine';
%  finelevel= '-fine';
   finelevel= '';
end

   ldpath='';
   if  exist('OCTAVE_VERSION') % FIXME
     islinux =1;
   elseif strfind(system_dependent('getos'),'Linux')
     islinux =1;
     s=version; ff= find(s=='.');
      if str2num(s(1:ff(2)-1))>=7
        %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks netgen    
          ldpath ='LD_LIBRARY_PATH=/usr/lib/Togl1.7:/opt/netgen/lib;';
      end      
   else
     islinux =0;
   end    

% Netgen executable filename
   if  islinux
      ng_name = '/opt/netgen/bin/netgen';
   else
      ng_name = 'ng';
   end

while( 1 )

   fid= fopen('ng.opt','w'); %create ng.opt file in local dir
   if fid==-1
      error(['Netgen requires writing files in the current directory(%s).', ...
             'Unfortunately, you don''t have permission'], pwd);
   end
   if ~isempty(msz_file)
%     fprintf(fid,'options.segmentsperedge 5\n'); % Another
%                                                   potentially useful parameter
%                                                   except netgen ignores it
      fprintf(fid,'options.meshsizefilename= %s\n',msz_file);
   end
   fclose(fid);

   status= system(sprintf( ...
        '%s %s %s -batchmode -geofile=%s  -meshfile=%s ', ...
         ldpath, ng_name, finelevel,geo_file,vol_file));
   if status==0; break; end

   if islinux
      disp(['It seems you are running Linux and netgen has not worked. ' ...
           'Check that it is installed and on the path. ' ...
           'Perhaps LD_LIBRARY_PATH needs to be set?' ...
           'For newer versions of netgen, you need to set the environment variable NETGENDIR.', ...
           'Please enter a new netgen file name' ]);
      ng_name = input('netgen file name (with path)? [or i=ignore, e=error] ','s');
      if strcmp(ng_name,'i'); break;end
      if strcmp(ng_name,'e'); error('user requested'),end;
   else
      fprintf([ ...
       'Netgen call failed. Is netgen installed and on the search path?\n' ...
       'If you are running under windows, I can attempt to create\n' ...
       'a batch file to access netgen.\n' ...
       'Please enter the directory in which to find netgen.exe.\n' ...
       'If you don''t have a copy, download it from' ...
       'http://sourceforge.net/projects/netgen-mesher/ \n\n' ...
       'Note that you *MUST* use names without spaces. Thus\n' ...
       'instead of C:/Program Files/ write C:/Progra~1/\n\n' ]);
      netgen_path = input('netgen_path? [or i=ignore, e=error] ','s');
      if strcmp(ng_name,'i'); break;end
      if strcmp(ng_name,'e'); error('user requested'),end;
      if exist( sprintf('%s/netgen.exe',netgen_path) , 'file' ) || ...
              exist( sprintf('%s/bin/netgen.exe',netgen_path) , 'file' )
         disp('Found netgen version 4.4 or higher');
         netgen_exe = netgen_path;
         if exist( sprintf('%s/bin/netgen.exe',netgen_path) , 'file' )
             netgen_exe = [netgen_path '/bin'];
         end

         
         fid= fopen('ng.bat','w');
         fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path); % REQ for ng <= 4.4
         fprintf(fid,'set TIX_LIBRARY=%s/lib/tix8.1\n', netgen_path); % REQ for ng <= 4.4
         fprintf(fid,'set NETGENDIR=%s\n', netgen_path); % REQ for ng >= 4.9
         fprintf(fid,'%s/netgen.exe %%*\n', netgen_exe);
         fclose(fid);
      elseif exist( sprintf('%s/ng431.exe',netgen_path) , 'file' ) 
         disp('Found netgen version 4.3.1');
         fid= fopen('ng.bat','w');
         fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
         fprintf(fid,'set TIX_LIBRARY=%s/lib/tcl8.2\n', netgen_path);
         fprintf(fid,'%s/ng431.exe %%*\n', netgen_path);
         fclose(fid);
      else
         warning(['cannot find a version of netgen that I know about\n' ...
                  'Install netgen 4.4 or 4.3.1 or check the path\n']);
      end
   end
end
