function print_convert(filename, varargin);
% PRINT_CONVERT: print figures and call imagemagick to trim them
% print_convert(filename, options, pagehwr)
%
%  filename: filename to print to file type is the extension
%  options:  options to send to imagemagick convert (if it is available)
%  pagehwr:  set page high/width to control shape of figures
%
%  Examples
%   print_convert('outname.png') % assumes resolution = 125dpi
%   print_convert('outname.png','-r150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)
%
% Compatibility with pre-3.7 features
%   print_convert('outname.png','-density 150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)

% More options
%  - options to the imwrite command
%  - options to the print command
%  - options - crop all spaces
%  - options - change page size

% (C) Andy Adler 2010. License: GPL v2 or v3.
% $Id$

if isstr(filename) && strcmp(filename,'UNIT_TEST'); do_unit_test; return; end

pp = parse_options(varargin{:});

tmpnam = [tempname,'.png'];

posn = get(gcf,'PaperPosition');
% I wish matlab gave us unwind protect - like octave does!
set(gcf,'PaperPosition',[posn(1:3), posn(3)*pp.pagehwr]);
print('-dpng',tmpnam);
set(gcf,'PaperPosition',[posn(1:4)]);

im = imread(tmpnam,'png');
delete(tmpnam);

im= crop_image(im,pp);
imwrite(im,filename);

function im = crop_image(im,pp)
   szim = size(im);
   bdr = mean(double(im(1,:,:)),2);

   isbdr = true(szim(1),szim(2));
   for i=1:szim(3);
     isbdr = isbdr & (im(:,:,i) == bdr(i));
   end

   horz = [true,all(isbdr,1),true];
   horzpt = find(diff(horz)) - 1; % first 'true'
   im(:,1:horzpt(1)    ,:)= [];
   im(:,horzpt(end):end,:)= [];

   vert = [true,all(isbdr,2)',true];
   vertpt = find(diff(vert)) - 1; % first 'true'
   im(1:vertpt(1)    ,:,:)= [];
   im(vertpt(end):end,:,:)= [];

function pp = parse_options(varargin)
   if nargin<=1; pp.options = '';  end
   if nargin<=2; pp.pagehwr = 6/8; end

function do_unit_test
   fid = fopen('print_convert_test.html','w');
   fprintf(fid,'<HTML><BODY>\n');
   for i=1:26;
     fprintf(fid,'<H1>%02d</H1><img src="pc%02d.png">\n',i,i);
   end
   fprintf(fid,'</BODY></HTML>\n');
   fclose(fid);

   im = mk_image( mk_common_model('b2c2',8), 1:256);
   clf; show_fem(im);
   print_convert pc01.png

   clf; subplot(221); show_fem(im); 
        subplot(224); show_slices(im); 
   print_convert pc02.png

   im = mk_image( ng_mk_cyl_models(1,[16,.5],.05),1);
   clf; 
   show_fem(im);
   print_convert pc03.png

   
