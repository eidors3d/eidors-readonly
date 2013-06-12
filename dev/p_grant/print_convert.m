function print_convert(filename, varargin)
% PRINT_CONVERT: print figures and trim them
% print_convert(filename, options, pagehwr)
%
%  filename: filename to print to file type is the extension
%  options:  options to send to imagemagick convert (if it is available)
%  pagehwr:  set page high/width to control shape of figures
%
%  Examples
%   print_convert('outname.png') % assumes resolution = 125dpi
%   print_convert('outname.png',opts);
%
%  Where options are:
%   opt.resolution = 150; % 150 dpi
%   opt.pagesize   = [width, height]; %whatever matlab's current units
%   opt.imwrite_opts = {'Quality',75}; % Jpeg quality
%   opt.supersampling_factor = 2; % default, moderate anti-aliasing
%
%  Other options are in the UNIT_TEST section of this file
%
% Compatibility with pre-3.7 features
%   print_convert('outname.png','-density 150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)
 
% (C) Andy Adler 2010-2013. License: GPL v2 or v3.
% $Id$
 
if ischar(filename) && strcmp(filename,'UNIT_TEST'); do_unit_test; return; end

pp = parse_options(varargin{:});

tmpnam = [tempname,'.png'];

old_ihc = get(gcf, 'InvertHardcopy');
old_col = get(gcf, 'Color');
set(gcf,'InvertHardCopy','off'); 
set(gcf,'Color','w');

set(gcf,'PaperPosition',pp.posn); % I wish matlab had unwind protect - like octave does!
print('-dpng ',pp.resolution,tmpnam);
set(gcf,'PaperPosition',pp.page);
 
set(gcf,'InvertHardCopy',old_ihc); % 
set(gcf,'Color',old_col);
 
im = imread(tmpnam,'png');
delete(tmpnam);

im = bitmap_downsize(im, pp.factor);
im = crop_image(im,pp);

if (strfind(filename,'.jpg') == (length(filename)-3) );
    imwrite(im,filename,'Quality', pp.jpeg_quality, pp.imwrite_opts{:});
else 
imwrite(im,filename,pp.imwrite_opts{:});
end 
 
function im = crop_image(im,pp)
   szim = size(im);
   bdr = squeeze(mean(double(im(1,:,:)),2));
 
   isbdr = true(szim(1),szim(2));
   for i=1:szim(3);
     isbdr = isbdr & (im(:,:,i) == bdr(i));
   end
 
   horz = [true,all(isbdr,1),true];
   horzpt = find(diff(horz)) - 1; % first 'true'
   im(:,horzpt(end)+1:end,:)= []; % remove higher first
   if pp.horz_cut >0;
      horz_e_pt = find(diff(horz)==-1) -1; horz_e_pt(1) = [];
      horz_s_pt = find(diff(horz)==+1)   ; horz_s_pt(end) = [];
      idx = find(horz_e_pt - horz_s_pt > pp.horz_cut);
      for i=fliplr(idx) % remove higher first
        im(:,horz_s_pt(i):horz_e_pt(i),:)= [];
      end
   end
   im(:,1:horzpt(1)    ,:)= [];
 
   vert = [true,all(isbdr,2)',true];
   vertpt = find(diff(vert)) - 1; % first 'true'
   im(vertpt(end)+1:end,:,:)= [];
   if pp.vert_cut >0;
      vert_e_pt = find(diff(vert)==-1) -1; vert_e_pt(1) = [];
      vert_s_pt = find(diff(vert)==+1)   ; vert_s_pt(end) = [];
      idx = find(vert_e_pt - vert_s_pt > pp.vert_cut);
      for i=fliplr(idx) % remove higher first
        im(vert_s_pt(i):vert_e_pt(i),:,:)= [];
      end
   end
   im(1:vertpt(1)    ,:,:)= [];
 
    
% factor = 1 if all plots only contain images, 2 otherwise
function f = default_factor
   f = 1; % default for images
   sp = get(gcf,'Children'); % subplots
   for i = 1:length(sp)
      obj = get(sp(i),'Children');
      tp  = get(obj,'Type');
      if ~all(strcmp(tp,'image'))
         f = 2;
         return;
      end
   end 
   
function pp = parse_options(varargin)
   pp.page = get(gcf,'PaperPosition');
   pp.posn = pp.page;
   pp.imwrite_opts = {};
   pp.horz_cut = 0;
   pp.vert_cut = 0;
   pp.factor = default_factor;
   pp.resolution = sprintf('-r%d',125 * pp.factor);
   pp.jpeg_quality = 80;
 
% Old options
   if nargin< 1; pp.options = '';  return; end
   if nargin>=2; pp.posn(4) = pp.posn(3)*varargin{2};  end
 
   opt = varargin{1};
   if ischar(opt)
      val =regexp(opt,'-density (\d+)','tokens');
      if length(val)>0;
         pp.resolution = sprintf('-r%d', str2double(val{1}{1}) * pp.factor);
      end
      val =regexp(opt,'-r(\d+)','tokens');
      if length(val)>0;
         pp.resolution = sprintf('-r%d', str2double(val{1}{1}) * pp.factor);
      end
   elseif isstruct(opt)
     if isfield(opt,'supersampling_factor')
        pp.factor = opt.supersampling_factor;
     end
     if isfield(opt,'resolution');
         pp.resolution = sprintf('-r%d', opt.resolution * pp.factor);
     else
         pp.resolution = sprintf('-r%d',125 * pp.factor);
     end
     if isfield(opt,'pagesize');
         pp.posn(3:4) = opt.pagesize;
     end
     % TODO, this code can copy from opt to pp
     if isfield(opt,'imwrite_opts');
        pp.imwrite_opts = opt.imwrite_opts;
        %if ~any(strcmp(pp.imwrite_opts,'quality') | strcmp(pp.imwrite_opts,'Quality') );
           %pp.imwrite_opts(end+1:end+2) = {'quality',pp.jpeg_quality};
        %end
     %else 
        %pp.imwrite_opts(end+1:end+2) = {'quality',pp.jpeg_quality};
     end
    
     if isfield(opt,'horz_cut');
         pp.horz_cut = opt.horz_cut;
     end
     if isfield(opt,'vert_cut');
         pp.vert_cut = opt.vert_cut;
     end
   else
      error('Can''t parse options');
   end
 
 
function do_unit_test
   fprintf('does unit test \n');
   fid = fopen('print_convert_test.html','w');
   fprintf(fid,'<HTML><BODY>\n');
   for i=1:26;
     switch i;
       case {9,10}; typ = 'jpg';
       otherwise; typ = 'png';
     end
     fprintf(fid,...
     '<H1>%02d</H1><table border=1><tr><td><img src="pc%02d.%s"></table>\n',i,i,typ);
   end
   fprintf(fid,'</BODY></HTML>\n');
   fclose(fid);
   eidors_msg('TO VIEW OUTPUT, OPEN FILE print_convert_test.html',1);
 
   im = mk_image( ng_mk_cyl_models(1,[16,.5],.05),1);
   clf; show_fem(im);
   print_convert pc01.png
 
   im = mk_image( mk_common_model('b2c2',8), 1:256);
   clf; show_fem(im);
   print_convert pc02.png
 
 
   clf; subplot(221); show_fem(im); 
        subplot(224); show_slices(im); 
   print_convert pc03.png
 
   print_convert pc04.png '-density 50'
   print_convert pc05.png '-r100'
 
   print_convert('pc06.png', '-r100' , 1/8)
 
   print_convert('pc06.png', '-r100' , 3/1)
 
   opt.resolution = 100;
   print_convert('pc07.png', opt);
   opt.pagesize = [8,2];
   print_convert('pc08.png', opt);
 
   print_convert('pc09a.jpg');
   
   opt.pagesize = [8,6];
   print_convert('pc09.jpg', opt);
  
   
   opt.imwrite_opts = {'Quality',20};
   print_convert('pc10.jpg', opt);
   
   opt.imwrite_opts = {'Quality',100};
   print_convert('pc10b.jpg', opt);
 
   opt.imwrite_opts = {};
   opt.horz_cut = 50;
   print_convert('pc11.png', opt);
 
   opt.vert_cut = 50;
   print_convert('pc12.png', opt);
 
   clf; subplot(331); show_fem(im); 
        subplot(335); show_slices(im); 
        subplot(339); show_slices(im); 
   opt.vert_cut = 30;
   print_convert('pc13.png', opt);
   
   im = mk_image( ng_mk_cyl_models(1,[16,.5],.05),1);
   clf; show_fem(im);
   clear opt; opt.supersampling_factor = 1;
   print_convert('pc14.png', opt);
 
   clear opt; opt.supersampling_factor = 2;
   print_convert('pc15.png', opt);
 
   clear opt; opt.supersampling_factor = 3;
   print_convert('pc16.png', opt);
   
   clear opt; opt.supersampling_factor = 4;
   print_convert('pc17.png', opt);
   
   clear opt; opt.supersampling_factor = 8;
   print_convert('pc18.png', opt);
 

