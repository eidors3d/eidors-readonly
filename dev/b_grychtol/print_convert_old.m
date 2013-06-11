function print_convert(filename, options, pagehwr)
% PRINT_CONVERT: print figures and call imagemagick to trim them
% print_convert(filename, options, pagehwr)
%
%  filename: filename to print to file type is the extension
%  options:  options to send to imagemagick convert (if it is available)
%  pagehwr:  set page high/width to control shape of figures
%
%  Examples
%   print_convert('outname.png') % assumes resolution = 125dpi
%   print_convert('outname.png','-density 150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)

% (C) Andy Adler 2010. License: GPL v2 or v3.
% $Id$

if nargin<=1; options = '';  end
if nargin<=2; pagehwr = 6/8; end

% Test if we have imagemagick convert
[ret1,ret2] =system('convert -version');
if ret1 == 1 || ~any(findstr(ret2, 'ImageMagick'))
   eidors_msg(['This function requires the ImageMagick convert program ' ...
               ' be installed and on the path.' ...
               ' On windows, we recommend you install it as part of ' ...
               ' cygwin (www.cygwin.com). Continuing without this sw.'],1);
   ext = find(filename=='.'); ext = ext(end);
   try;  ext = filename(ext+1:end);
   catch; error('unable to get filename extension');
   end
   if strcmp(ext,'jpg'); ext='jpeg';  end
   if strcmp(ext,'ps');  ext='psc2';  end
   if strcmp(ext,'eps'); ext='epsc2'; end

   posn = get(gcf,'PaperPosition');
   set(gcf,'PaperPosition',[posn(1:3), posn(3)*pagehwr]);
   print(['-d',ext],filename);
   set(gcf,'PaperPosition',[posn(1:4)]);
   return
end

tmpnam = [tempname,'.eps'];

posn = get(gcf,'PaperPosition');
% I wish matlab gave us unwind protect - like octave does!
set(gcf,'PaperPosition',[posn(1:3), posn(3)*pagehwr]);
print('-depsc2',tmpnam);
set(gcf,'PaperPosition',[posn(1:4)]);

% Fix cygwin bugs
slash = find(tmpnam=='\');
if any(slash)
   tmpnam(slash) = '/'; % Works for all oses
%  tmpnam = ['/cygdrive/', tmpnam([1,3:end])]; - find a way to generalize
end


% this isn't working for old versions of vnc
% Need colorspace RGB, otherwise IE can't read it.
% 2012/12/6: RGB messes up the colors, sRGB does not
cmd = sprintf('convert -colorspace sRGB -density 125 %s -trim  %s %s',  options, tmpnam, filename);

% Code to try to use gs to get around imagemagick bugs
if 0
options = regexprep(options, '\-density (\d+)','-r$1');
cmd = sprintf('gs -r125 %s -dEPSCrop -sDEVICE=png16m -sOutputFile=%s -dBATCH -dNOPAUSE %s', ...
     options, filename, tmpnam) 
end

flag = system_cmd(cmd);
if flag~=0
   warning('Failed to call imagemagick to crop image?');
end
