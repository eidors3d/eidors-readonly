function print_convert(filename, varargin)
%PRINT_CONVERT: print figures with anti-aliasing and trim them
%  PRINT_CONVERT(FILENAME,OPT) prints current figure to FILENAME, which 
%  must include extenstion.
%
%  OPT is either a string specifying dpi in the format '-r150' or a struct
%  with (some of) these fields:
%   opt.resolution   = 150;              % 150 dpi (default: 125)
%   opt.pagesize     = [width, height];  % in current PageUnits or 'auto'
%                                        % to match the size on the screen
%                                        % (default: [8, 6] inches)
%   opt.jpeg_quality = 85;               % jpeg quality (default: 85)
%   opt.imwrite_opts = {'BitDepth',8};   % options to IMWRITE (default: '')
%   opt.horz_cut     = 20;               % remove horizontal gaps larger
%                                        % than 50 pixels (default: 20)
%   opt.horz_space   = 10;               % replace the removed horizontal
%                                        % gaps with 10 px of background
%                                        % (default: 10)
%   opt.vert_cut     = 20;               % remove vertical gaps larger
%                                        % than 50 pixels (default: 20)
%   opt.vert_space   = 10;               % replace the removed vertical
%                                        % gaps with 10 px of background
%                                        % (default: 10)
%   opt.supersampling_factor = 2;        % anti-aliasing (default: 1 for
%                                        % images, 2 for graphs). Higher is
%                                        % smoother.
%   opt.crop_slack   = [top,bot,left,right] %don't crop right to boundary
%   opt.figno        = number            % print figure, otherwise gca
%
%  Note that opt.imwrite_opts takes precedence over opt.jpeq_quality.
%
%  print_convert has pre-version 3.7 options, which are deprecated
%
%  Examples
%   print_convert('outname.png')         % uses default options
%   print_convert outname.png
%   print_convert outname.png -r150       % set dpi=150
%   print_convert('outname.png',opts);   % use options
%
% See also IMWRITE
 
% (C) Andy Adler and Bartlomiej Grychtol 2010-2021. License: GPL v2 or v3.
% $Id$
%
%  Compatibility with pre-3.7 features:
%  ------------------------------------
%  PRINT_CONVERT(FILENAME, OPTIONS, PAGEHWR)
%
%  FILENAME : filename to print to file type is the extension
%  OPTIONS  : specify dpi as '-density N'
%  PAGEHWR  : set page hight/width to control shape of figures
%  
%   print_convert('outname.png','-density 150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)
 
if ischar(filename) && strcmp(filename,'UNIT_TEST'); do_unit_test; return; end

pp = parse_options(filename, varargin{:});


% change properties 
pp = set_prop(pp,'InvertHardCopy','off');
pp = set_prop(pp,'Color','w');
pp = set_pagesize(pp);

% revert changed properties at exit 
cleanupObj = onCleanup(@() set(pp.figno, pp.old{:})); % current state of pp.old only!

% print to array
im = print(pp.figno,'-RGBImage',pp.resolution);

im = bitmap_downsize(im, pp.factor);
im = crop_image(im,pp);
try
   imwrite(im,pp.filename,pp.imwrite_opts{:});
catch e
   eidors_msg(['Call to IMWRITE failed.'...
               'Probably opt.imwrite_opts is incorrect for %s files.'],...
               upper(pp.fmt), 1);
   disp('opt.imwrite_opts:');
   disp(pp.imwrite_opts);
   rethrow(e);
end

% place file in clipboard
if pp.clip
   copy2clipboard(pp)
end

function pp = set_pagesize(pp)
if ~isfield(pp, 'pagesize') % no user definition
    pos = [0.25 0.25 8 6]; % Matlab's default before R2016a
    if isfield(pp, 'aspect_ratio') % old interface
        pos(4) = pos(3) * pp.aspect_ratio;
    end
    pp = set_prop(pp,   'PaperUnits','inches',...
                        'PaperPositionMode','manual', ...
                        'PaperPosition', pos);
elseif ischar(pp.pagesize)
    switch pp.pagesize
        case 'auto'
            pp = set_prop(pp, 'PaperPositionMode', 'auto');
        otherwise
            warning('Page size %s not understood', pp.pagesize);
    end
else
    pos = [0 0 pp.pagesize];
    pp = set_prop(pp,   'PaperPositionMode','manual', ...
                        'PaperPosition', pos);
end    

function pp = set_prop(pp, varargin)
    prop = {};
    for i = 1:2:numel(varargin)
        prop(1, end+1) = varargin(i);
        prop{2, end  } = get(pp.figno, varargin{i});
    end
    pp.old = [fliplr(prop), pp.old]; % store in reverse order to revert correctly
    set(pp.figno, varargin{:})
        

function copy2clipboard(pp)
   if isunix()
      cmd = sprintf(['xclip -selection ' ... 
          'clipboard -t image/%s -i %s'], ...
            pp.fmt, pp.filename);
      [status] = system(cmd);
      if status ~=0
          warning 'calling xclip to clipboard failed';
      end
   else
      cmd = 'Powershell -command Set-Clipboard -Path ';
      status = system([cmd pp.filename]);
      if status ~=0
        warning 'Powershell Set-Clipboard failed';
      end
   end

function im = crop_image(im,pp)
   tu = pp.crop_slack(1);
   bu = pp.crop_slack(2) + 1; %remove starting one more
   lu = pp.crop_slack(3);
   ru = pp.crop_slack(4) + 1;

   szim = size(im);
   bdr = squeeze(mean(double(im(1,:,:)),2));
 
   isbdr = true(szim(1),szim(2));
   for i=1:szim(3);
     isbdr = isbdr & (im(:,:,i) == bdr(i));
   end
 
   horz = [true,all(isbdr,1),true];
   horzpt = find(diff(horz)) - 1; % first 'true'
   if isempty(horzpt)
      eidors_msg('Image is blank. Cropping aborted.',1);
      return
   end
   im(:,horzpt(end)+ru:end,:)= []; % remove higher first
   if pp.horz_cut >0;
      horz_e_pt = find(diff(horz)==-1) -1; horz_e_pt(1) = [];
      horz_s_pt = find(diff(horz)==+1)   ; horz_s_pt(end) = [];
      idx = find(horz_e_pt - horz_s_pt > pp.horz_cut);
      for i=fliplr(idx) % remove higher first
        im(:,horz_s_pt(i)+pp.horz_space:horz_e_pt(i),:)= [];
      end
   end
   im(:,1:horzpt(1)-lu ,:)= [];
 
   vert = [true,all(isbdr,2)',true];
   vertpt = find(diff(vert)) - 1; % first 'true'
   im(vertpt(end)+bu:end,:,:)= [];
   if pp.vert_cut >0;
      vert_e_pt = find(diff(vert)==-1) -1; vert_e_pt(1) = [];
      vert_s_pt = find(diff(vert)==+1)   ; vert_s_pt(end) = [];
      idx = find(vert_e_pt - vert_s_pt > pp.vert_cut);
      for i=fliplr(idx) % remove higher first
        im(vert_s_pt(i)+pp.vert_space:vert_e_pt(i),:,:)= [];
      end
   end
   im(1:vertpt(1)-tu ,:,:)= [];
 
    
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
   
function fmt = parse_format(filename)   
    ext = lower(regexp(filename,'(?<=\.).+$','match'));
    if isempty(ext); error('no filename extension detected (%s)',filename); end
    switch ext{1}
       case {'jpg', 'jpeg'}
          fmt = 'jpg';
       case {'j2c', 'j2k', 'jp2'}
          fmt = 'jp2';
       case {'tif','tiff'}
          fmt = 'tif';
       otherwise
          fmt = ext{1};
    end
       
function pp = parse_options(filename,varargin)
   
   pp.fmt = parse_format(filename);
   if strcmp(filename, ['{clipboard}.',pp.fmt])
      pp.clip = true;
      pp.filename = [tempname(),'.',pp.fmt]; 
   else
      pp.clip = false;
      pp.filename = filename;
   end

   pp.jpeg_quality = 85; % default jpeg quality
   pp.imwrite_opts = {}; % initial
   pp.horz_cut = 20;
   pp.horz_space = 10;
   pp.vert_cut = 20;
   pp.vert_space = 10;
   pp.factor = default_factor;
   pp.resolution = sprintf('-r%d',125 * pp.factor);
   pp.crop_slack = [0,0,0,0];
   if exist('OCTAVE_VERSION')
      pp.figno = gcf;
   else
      pp.figno = get(gcf,'Number');
   end
   pp.old = {};
   
 
% Old options
   if nargin == 1   
      return; 
   end
   if nargin>=3
      %print_convert('outname.png','-density 150', 0.5)
      pp.aspect_ratio = varargin{2};  
   end
 
   opt = varargin{1};
   if ischar(opt)
      val =regexp(opt,'-density (\d+)','tokens');
      if ~isempty(val);
         pp.resolution = sprintf('-r%d', str2double(val{1}{1}) * pp.factor);
      end
      val =regexp(opt,'-r(\d+)','tokens');
      if ~isempty(val);
         pp.resolution = sprintf('-r%d', str2double(val{1}{1}) * pp.factor);
      end
   elseif isstruct(opt)
     if isfield(opt,'figno');
        pp.figno = opt.figno;
     end
     if isfield(opt,'supersampling_factor')
        pp.factor = opt.supersampling_factor;
     end
     if isfield(opt,'resolution');
         pp.resolution = sprintf('-r%d', opt.resolution * pp.factor);
     else
         pp.resolution = sprintf('-r%d',125 * pp.factor);
     end
     if isfield(opt,'pagesize');
         pp.pagesize = opt.pagesize;
     end
     % TODO, this code can copy from opt to pp
     if isfield(opt,'jpeg_quality')
        pp.jpeg_quality = opt.jpeg_quality;
     end
     if strcmp(pp.fmt,'jpg')
        pp.imwrite_opts = {'quality',pp.jpeg_quality}; 
     end
     if isfield(opt,'imwrite_opts');
        pp.imwrite_opts = opt.imwrite_opts;
        if strcmp(pp.fmt,'jpg') && ~any(strcmpi(pp.imwrite_opts,'quality'))
           pp.imwrite_opts(end+1:end+2) = {'quality',pp.jpeg_quality};
        end
     end
    
     if isfield(opt,'horz_cut');
         pp.horz_cut = opt.horz_cut;
     end
     if isfield(opt,'vert_cut');
         pp.vert_cut = opt.vert_cut;
     end
     if isfield(opt,'vert_space');
        if opt.vert_space >= pp.vert_cut;
           warning('Option vert_space must be smaller than vert_cut. Ingoring');
        else
           pp.vert_space = opt.vert_space;
        end
     end
     if isfield(opt,'horz_space');
        if opt.horz_space >= pp.horz_cut;
           warning('Option vert_space must be smaller than vert_cut. Ingoring');
        else
           pp.horz_space = opt.horz_space;
        end
     end
     if isfield(opt,'crop_slack');
        pp.crop_slack = opt.crop_slack;
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
 

