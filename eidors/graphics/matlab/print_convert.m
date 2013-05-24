function print_convert(filename, varargin)
% PRINT_CONVERT - print figures and trim them (export_fig wrapper)
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
%
%  Other options are in the UNIT_TEST section of this file
%
% Compatibility with pre-3.7 features
%   print_convert('outname.png','-density 150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)
%
% NOTE: pagehwr has no effect at the moment


% (C) Andy Adler 2010-2013. License: GPL v2 or v3.
% $Id$



pp = parse_options(varargin{:});


old_col = get(gcf, 'Color');
set(gcf,'InvertHardCopy','off'); % 
set(gcf,'Color','w');

set(gcf,'PaperPosition',pp.posn); % I wish matlab had unwind protect - like octave does!
opt = sprintf('-r%s', pp.resolution);
export_fig(filename,opt)

set(gcf,'PaperPosition',pp.page);
set(gcf,'Color',old_col);


function pp = parse_options(varargin)
   pp.resolution = '125';
   pp.page = get(gcf,'PaperPosition');
   pp.posn = pp.page;

% Old options
   if nargin< 1; pp.options = '';  return; end
   if nargin>=2; pp.posn(4) = pp.posn(3)*varargin{2};  end

   opt = varargin{1};
   if isstr(opt)
      val =regexp(opt,'-density (\d+)','tokens');
      if length(val)>0;
         pp.resolution = ['-r',val{1}{1}];
      end
      val =regexp(opt,'-r(\d+)','tokens');
      if length(val)>0;
         pp.resolution = ['-r',val{1}{1}];
      end
   elseif isstruct(opt)
     if isfield(opt,'resolution');
         pp.resolution = sprintf('-r%d',opt.resolution);
     end
     if isfield(opt,'pagesize');
         pp.posn(3:4) = opt.pagesize;
     end
   else
      error('Can''t parse options');
   end