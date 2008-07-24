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
% $Id: call_netgen.m,v 1.9 2008-06-11 19:19:29 aadler Exp $
% (C) 2006 Andy Adler. Licensed under GPL V2

if nargin<3
   msz_file= '';
end

if nargin<4
%  finelevel= '-veryfine';
%  finelevel= '-fine';
   finelevel= '';
end

while( 1 )
   ldpath='';
   if  strfind(system_dependent('getos'),'Linux')
     islinux =1;
     s=version; ff= find(s=='.');
      if str2num(s(1:ff(2)-1))>=7
        %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks netgen    
          ldpath ='LD_LIBRARY_PATH=;';
      end      
   else
     islinux =0;
   end    

   fid= fopen('ng.opt','w'); %create ng.opt file in local dir
   if ~isempty(msz_file)
%     fprintf(fid,'options.segmentsperedge 5\n'); % Another
%                                                   potentially useful parameter
%                                                   except netgen ignores it
      fprintf(fid,'options.meshsizefilename  %s\n',msz_file);
   end
   fclose(fid);

   status= system(sprintf( ...
        '%s ng %s -batchmode -geofile=%s  -meshfile=%s ', ...
         ldpath,finelevel,geo_file,vol_file));
   if status==0; break; end

   if islinux
       disp('It seems you are running Linux and netgen has not worked. Check that it is installed and on the path. Perhaps LD_LIBRARY_PATH needs to be set?');
   else

   fprintf([ ...
    'Netgen call failed. Is netgen installed and on the search path?\n' ...
    'If you are running under windows, I can attempt to create\n' ...
    'a batch file to access netgen.\n' ...
    'Please enter the directory in which to find netgen.\n' ...
    'If you don''t have a copy, download it from' ...
    'http://www.hpfem.jku.at/netgen/\n\n' ...
    'Note that you *MUST* use names without spaces. Thus\n' ...
    'instead of C:/Program Files/ write C:/Progra~1/\n\n' ...
    ]);
   netgen_path = input('netgen_path? ','s');
   if exist( sprintf('%s/netgen.exe',netgen_path) , 'file' ) 
      disp('Found netgen version 4.4');

      fid= fopen('ng.bat','w');
      fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
      fprintf(fid,'set TIX_LIBRARY=%s/lib/tix8.1\n', netgen_path);
      fprintf(fid,'%s/netgen.exe %%*\n', netgen_path);
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
