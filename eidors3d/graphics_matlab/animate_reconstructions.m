function fname_out= animate_reconstructions(fname, imgs);
% animate_reconstructions(fname, imgs);
% animate a sequence of reconstructed images
%
% PARAMETER:  fname
%     filename to save to (extension is added)
% PARAMETER:  imgs
%     is array of eidors images
%
% OUTPUT: fname_out
%     Name of animated file written to.
%     An animated window will not pop up if output requested

mk_movie2(fname,imgs)

if nargout>0
   fname_out= [fname, '.gif'];
else
   fid= fopen([fname ,'.html'],'w');
   fprintf(fid,'<HTML><HEAD><TITLE>%s</TITLE><BODY>\n',fname);
   fprintf(fid,'<img src="%s.gif" width="256"></BODY></HTML>',fname);
   fclose(fid);
   if  strfind(system_dependent('getos'),'Windows')
      system(sprintf('explorer "%s.html"',fname));
   else % we hope this is here - under linux etc
      system(sprintf('mozilla "./%s.html"',fname));
   end
end

% create avi movie fname
% imds is array of eidors images
%
% This results in really poor images with lots
%  of compression artefacts. NO really good reason to use
function mk_movie(fname, imgs, clim, ref_lev)
   fig=figure;
   set(fig,'DoubleBuffer','on');
   mov = avifile( [fname ,'.avi'] );%, 'Compression', 'RLE' );

   for i=1:length(imgs)
     show_slices(imgs(i),1,clim,ref_lev);
     F = getframe(gca);
     mov = addframe(mov,F);
   end
   mov = close(mov);
   close(fig);

% create gif movie fname
% imgs is array of eidors images
%
% This requires imagemagick convert program.
function mk_movie2(fname, imgs);
   calc_colours('mapped_colour', 127);
   dirname= 'tmp_mk_movie2';
   rm_rf( dirname );
   mkdir( dirname );

   r_img= calc_slices(imgs);
   c_img = calc_colours( r_img, imgs);
   out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
   cmap= colormap;

   for i=1:size(out_img,3)
     imwrite(out_img(:,:,i),cmap, ...
            sprintf('%s/img%05d.png',dirname, i), 'png');
   end

   ld_lib_path= sys_dep;

   if 0 % enumerate each file
      files= dir(sprintf('%s/img*.png', dirname ));
      % font selected is a windows font - how to make os-neutral?
      for ff= files(:)'
         fn= [dirname,'/',ff.name];
         fno= ff.name(4:8);
         retval= system(sprintf( ...
          '%s convert -font 6x8 -draw "text 0,10 ''%s''" %s %s', ...
          ld_lib_path, fno, fn, fn ));
      end
   end
      
   retval= system(sprintf( ...
       '%s convert -delay 5 %s/img*.png -loop 0 %s.gif', ...
       ld_lib_path, dirname, fname ));
   if retval~=0
       error('please ensure the imagemagick convert program is in your path. Under windows the easist is to download from www.imagemagick.org/script/binary-releases.php');
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
     vv=version; 
     ff=find(vv == ' '); 
     if length(ff)>0;vv=vv(1:ff(1)-1);end
     ff=find(vv == '.');
     if length(ff)>1;vv=vv(1:ff(2)-1);end     
     if str2num(vv)>=7
        %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks external progs
          ld_lib_path='LD_LIBRARY_PATH=;';
      end      
   end    
