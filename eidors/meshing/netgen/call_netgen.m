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

if  exist('OCTAVE_VERSION') % FIXME
   islinux =1;
elseif ~strncmp(computer,'PC',2) % Don't know if we have isunix
   islinux =1;
else
   islinux =0;
end

% Netgen executable filename
cache_path = eidors_cache('cache_path');
if  islinux
   ng_name = '';
else
   ng_name = [cache_path,'/ng'];
end
 
while( 1 )
   
   fid= fopen('ng.opt','a'); %create ng.opt file in local dir
   if fid==-1
      error(['Netgen requires writing files in the current directory(%s). ', ...
         'Unfortunately, you don''t have permission. ' ...
         'Your options are: 1) change your working directory to one in which you have write permission, or ' ...
         '2) change the permissions on the current working directory.'], pwd);
   end
   if ~isempty(msz_file)
      %     fprintf(fid,'options.segmentsperedge 5\n'); % Another
      %                                                   potentially useful parameter
      %                                                   except netgen ignores it
      fprintf(fid,'options.meshsizefilename= %s\n',msz_file);
   end
   fclose(fid);

   if strncmp(computer,'PC',2)
      % on Linux, Netgen runs in the current directory
      % enforce this behavioud in Windows
      oldpath = getenv('NETGEN_USER_DIR');
      setenv('NETGEN_USER_DIR', cd);
   end

   if eidors_debug('query','call_netgen')
      sys_cmd = sprintf('"%s" %s  -geofile=%s  -meshfile=%s ', ...
         ng_name, finelevel,geo_file,vol_file);
   else
      sys_cmd = sprintf('"%s" %s -batchmode -geofile=%s  -meshfile=%s ', ...
         ng_name, finelevel,geo_file,vol_file);
   end
   status= system_cmd( sys_cmd );

   if status==0; break; end
   try
      if islinux
         fprintf(['It seems you are running Linux and netgen has not worked. ' ...
            'Check that it is installed and on the path. ' ...
            'For newer versions of netgen, you need to set the environment' ...
            ' variable NETGENDIR. This can be set by starting matlab as follows:\n', ...
            '   NETGENDIR=/path/to/netgen PATH=/path/to/netgen:$PATH matlab\n', ...
            'Please enter a new netgen file name\n' ]);
         ng_name = input('netgen file name (with path)? [or i=ignore, e=error] ','s');
         if strcmp(ng_name,'i'); break;end
         if strcmp(ng_name,'e'); error('user requested'),end;
      else
         fprintf([ ...
            'Netgen call failed. Is netgen installed and on the search path?\n' ...
            'If you are running under windows, I can attempt to create\n' ...
            'a batch file to access netgen.\n' ...
            'Please enter the directory in which to find netgen.exe.\n' ...
            'A typical path is "C:\\Program Files (x86)\\Netgen-5.0_x64\\bin"\n' ...
            'If you don''t have a copy, download it from' ...
            'http://sourceforge.net/projects/netgen-mesher/ \n\n']);
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
            
            
            fid= fopen([cache_path, '/ng.bat'],'w');
            if fid<0; error('Unable to write to %s',cache_path); end
            fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path); % REQ for ng <= 4.4
            fprintf(fid,'set TIX_LIBRARY=%s/lib/tix8.1\n', netgen_path); % REQ for ng <= 4.4
            fprintf(fid,'set NETGENDIR=%s\n', netgen_path); % REQ for ng >= 4.9
            fprintf(fid,'"%s/netgen.exe" %%*\n', netgen_exe);
            fclose(fid);
         elseif exist( sprintf('%s/ng431.exe',netgen_path) , 'file' )
            disp('Found netgen version 4.3.1');
            fid= fopen([cache_path, '/ng.bat'],'w');
            if fid<0; error('Unable to write to %s',cache_path); end
            fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
            fprintf(fid,'set TIX_LIBRARY=%s/lib/tcl8.2\n', netgen_path);
            fprintf(fid,'"%s/ng431.exe" %%*\n', netgen_path);
            fclose(fid);
         else
            warning(['cannot find a version of netgen that I know about\n' ...
               'Install netgen or check the path\n']);
         end
      end
   catch e
      if strncmp(computer,'PC',2)
         % restore Netgen settings on Windows
         setenv('NETGEN_USER_DIR', oldpath);
      end
      rethrow(e)
   end
end
