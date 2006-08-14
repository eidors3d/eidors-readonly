function animate_reconstructions(fname, imgs);
% animate_reconstructions(fname, imgs);
% animate a sequence of reconstructed images
% fname  is filename to save to (extension is added)
% imgs   is array of eidors images

mk_movie2(fname,imgs);
if  strfind(system_dependent('getos'),'Windows')
   fid= fopen([fname ,'.html'],'w');
   fprintf(fid,'<HTML><HEAD><TITLE>%s</TITLE><BODY>\n',fname);
   fprintf(fid,'<img src="%s.gif" width="256"></BODY></HTML>',fname);
   fclose(fid);
   system(sprintf('explorer "%s.html"',fname));
end

% create avi movie fname
% imds is array of eidors images
%
% This results in really poor images with lots
%  of compression artefacts. NO really good reason to use
function mk_movie(fname, imgs)
   fig=figure;
   set(fig,'DoubleBuffer','on');
   mov = avifile( [fname ,'.avi'] );%, 'Compression', 'RLE' );

   for i=1:length(imgs)
     show_slices(imgs(i));
     F = getframe(gca);
     mov = addframe(mov,F);
   end
   mov = close(mov);
   close(fig);

% create gif movie fname
% imgs is array of eidors images
%
% This requires imagemagick convert program.
function mk_movie2(fname, imgs)
   calc_colours('mapped_colour', 127);
   dirname= 'tmp_mk_movie2';
   rm_rf( dirname );
   mkdir( dirname );

   r_img= show_slices(imgs);
   c_img = calc_colours( r_img);
   out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
   cmap= colormap;

   for i=1:length(imgs)
     imwrite(out_img(:,:,i),cmap, ...
            sprintf('%s/img%05d.png',dirname, i), 'png');
   end

   ld_lib_path= sys_dep;

   retval= system(sprintf( ...
       '%s convert -delay 25 %s/img*.png -loop 0 %s.gif', ...
       ld_lib_path, dirname, fname ));
   if retval~=0
       error('please ensure the imagemagick convert program is in your path');
   end
   rm_rf(dirname);
   fprintf('file %s.gif created (in current directory)\n',fname);

function rm_rf(dirname)
   if isdir(dirname)==0
       return
   end

   if isunix
       system(['rm -rf "',dirname,'"']);
   else
       system(['rmdir /s /q "',dirname,'"']);
   end

% work around stupid matlab bugs
function ld_lib_path= sys_dep;
   ld_lib_path='';
   if  strfind(system_dependent('getos'),'Linux')
     s=ver;
      if str2num(s.Version)>=7
        %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks external progs
          ld_lib_path='LD_LIBRARY_PATH=;';
      end      
   end    
