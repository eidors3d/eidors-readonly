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

% (C) Andy Adler 2010. License: GPL v2 or v3. $Id$

if nargin<=1; options = '';  end
if nargin<=2; pagehwr = 6/8; end

% Test if we have imagemagick convert
[ret1,ret2] =system('convert -version');
if ret1 == 1 || ~any(findstr(ret2, 'ImageMagick'))
   eidors_msg(['This function requires the ImageMagick convert program ' ...
               ' be installed and on the path.' ...
               ' On windows, we recommend you install it as part of ' ...
               ' cygwin (www.cygwin.com)']);
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

ld = ''; % OVERRIDE STUPID MATLAB LD_LIBRARY_PATH
if isunix && ~exist('OCTAVE_VERSION','var');
   ld= 'LD_LIBRARY_PATH=""';
end

% this isn't working for old versions of vnc
cmd = sprintf('%s convert -density 125 %s -trim  %s %s', ld, options, tmpnam, filename);

% Code to try to use gs to get around imagemagick bugs
if 0
options = regexprep(options, '\-density (\d+)','-r$1');
cmd = sprintf('%s gs -r125 %s -dEPSCrop -sDEVICE=png16m -sOutputFile=%s -dBATCH -dNOPAUSE %s', ...
    ld, options, filename, tmpnam) 
end

flag = system(cmd);
if flag~=0
   warning('Failed to call imagemagick to crop image?');
end
