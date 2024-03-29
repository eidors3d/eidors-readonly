function fname_out= animate_reconstructions(fname, imgs);
% animate_reconstructions(fname, imgs);
% animate a sequence of reconstructed images
%
% PARAMETER:  fname
%     filename to save to (extension is added)
% PARAMETER:  imgs
%     is array of eidors images
%
% if imgs.animate_reconstructions.show_times = 1
%   then a timescale is shown on the bottom
% if imgs.animate_reconstructions.make_avi = 1
%   then use ffmpeg to write an avi
% if imgs.animate_reconstructions.frame_rate = 5 (images / sec)
%   default value is 5
%
% OUTPUT: fname_out
%     Name of animated file written to.
%     An animated window will not pop up if output requested

% (C) 2006 Andy Adler. Licensed under GPL version 2 or 3
% $Id$

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
      system(sprintf('firefox "./%s.html" &',fname));
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

   try
     show_times = imgs.animate_reconstructions.show_times;
   catch 
     show_times = 0;
   end

   try
     frame_rate = imgs.animate_reconstructions.frame_rate;
   catch 
     frame_rate = 5;
   end

   try 
     make_avi = imgs.animate_reconstructions.make_avi;
   catch 
     make_avi = 0;
   end

   if make_avi == 1; itp = 'jpg';
   else              itp = 'png';
   end

   r_img= calc_slices(imgs);
   c_img = calc_colours( r_img, imgs);
   out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
   cmap= colormap;

   [len_vi, len_hi, len_oi] = size(out_img);

   for i=1:len_oi
     this_img  = out_img(:,:,i);
     if show_times % add scrollbar on bottom
        add_bar = mk_add_bar( (i-1)/(len_oi-1), len_hi );
        this_img= [this_img; add_bar];
     end
     this_name = sprintf('%s/img%06d.%s',dirname, i, itp);
     imwrite(this_img, cmap, this_name, itp);
   end

   if 0 % enumerate each file
      files= dir(sprintf('%s/img*.%s', dirname,itp ));
      % font selected is a windows font - how to make os-neutral?
      for ff= files(:)'
         fn= [dirname,'/',ff.name];
         fno= ff.name(4:8);
         retval= system_cmd(sprintf( ...
          'convert -font 6x8 -draw "text 0,10 ''%s''" %s %s', ...
          fno, fn, fn ));
      end

   end
      
   if make_avi == 0
      retval= system_cmd(sprintf( ...
          'convert -delay %d %s/img*.png -loop 0 %s.gif', ...
          frame_rate, dirname, fname ));
   else
      retval= system_cmd(sprintf( ...
          'ffmpeg -r 25 -i %s/img*.jpg -vcodec msmpeg4v2 -y -an -qscale 2 %s.avi', ...
          dirname, fname ));
   end

   if retval~=0
       error('please ensure the imagemagick convert program is in your path. Under windows the easist is to download from www.imagemagick.org/script/binary-releases.php');
   end
   rm_rf(dirname);
   if make_avi == 0
   fprintf('file %s.gif created (in current directory)\n',fname);
   else
   fprintf('file %s.avi created (in current directory)\n',fname);
   end

function rm_rf(dirname)
   if isdir(dirname)==0
       return
   end

   if isunix
       system(['rm -rf "',dirname,'"']);
   else
       system(['rmdir /s /q "',dirname,'"']);
   end

function add_bar = mk_add_bar(frac, len) 
   sz_bar = 3/len;
   ind_val = 90;
   xax= linspace(0,1,len) - frac;
   yax= 1-abs( xax/sz_bar);
   yax= yax.*(yax>0);

   add_bar = round(ind_val * yax);
