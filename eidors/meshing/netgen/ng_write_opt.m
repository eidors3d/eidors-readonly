function opt = ng_write_opt(varargin)
%NG_WRITE_OPT Write an ng.opt file in current directory
%  NG_WRITE_OPT, without inputs, creates a default ng.opt
%
%  NG_WRITE_OPT(OPTSTRUCT) creates uses options as specified in OPTSTRUCT
%
%  NG_WRITE_OPT('sec.option',val,...) offers an alternative interface to
%  modify specific options
%
%  NG_WRITE_OPT(OPTSTR,...) where OPTSTR is one of 'very_coarse', 'coarse',
%  'moderate', 'fine', and 'very_fine', based on the corresponding defaults
%  in Netgen 5.0. Any additional inputs modify selected fields of that
%  default. NG_WRITE_OPT('moderate',...) is equivalent to
%  NG_WRITE_OPT(...).
%
%  OPTSTRUCT = NG_WRITE_OPT(...) returns a struct with the options.
%  No file is written.
%
%  NG_WRITE_OPT(PATH) copies the file specified in PATH as ng.opt to the
%  current directory.
%
% 
%  Note: some options only take effect when meshoptions.fineness is 6
%  (custom).
%
%  BEWARE: NG_WRITE_OPT will overwrite any existing ng.opt file in the current
%  directory. 
% 
%  Example:
%  ng_write_opt('meshoptions.fineness',6,'options.minmeshsize',20);
%  call_netgen(...)
%  delete('ng.opt'); % clean up
%
% To specify a volume for refinement, specify
%  ng_write_opt('MSZPOINTS', [list of x,y,z,maxh])
%   or
%  ng_write_opt('MSZBRICK', [xmin, xmax, ymin, ymax, zmin, zmax, maxh])
%
%  See also CALL_NETGEN

% (C) 2012 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

%TODO:
% 1. Check which options require meshoptions.fineness = 6 and enforce it

% if input is 'UNIT_TEST', run tests
if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'UNIT_TEST') 
   do_unit_test; return; end

if nargin == 1 && ischar(varargin{1}) &&  exist(varargin{1},'file') % a path
   copyfile(varargin{1},'ng.opt');
   return
end

nargs = nargin;

str = {};
% modify as per user input
if nargin >= 1 && ischar(varargin{1})
   str = varargin(1);
   varargin(1) = [];
   nargs = nargs - 1;
end

% get default options
opt = default_opt(str{1});

if nargs == 0
    usr = struct;
end
% modify as per user input
if nargs == 1 && isstruct(varargin{1})
   usr = varargin{1};
end
if nargs > 1 % string value pairs
   usr = process_input(varargin{:}); 
end
opt = copy_field(opt,usr);

if nargout == 1 % do not write a file if output was requested
   return
end

% write the file
write_ng_file(opt);

function opt = process_input(varargin)
if mod(nargin,2)~=0
   error('EIDORS:WrongInput','Even number of inputs expected');
end
for i = 1:nargin/2
   idx = (i-1)*2 +1;
   val = varargin{idx+1};
   switch varargin{idx}
%  ng_write_opt('MSZPOINTS', [list of x,y,z,maxh])
   case 'MSZPOINTS'
    %  ng_write_opt('MSZPOINTS', [list of x,y,z,maxh])
       tmpname = write_tmp_mszfile( val );
       opt.options.meshsizefilename = tmpname;
   case 'MSZBRICK'
    %  ng_write_opt('MSZBRICK', [xmin, xmax, ymin, ymax, zmin, zmax, maxh])
       maxh = val(7);
       xpts= floor(abs(diff(val(1:2))/maxh))+1;
       ypts= floor(abs(diff(val(3:4))/maxh))+1;
       zpts= floor(abs(diff(val(5:6))/maxh))+1;
       xsp= linspace(val(1),val(2), xpts);
       ysp= linspace(val(3),val(4), ypts);
       zsp= linspace(val(5),val(6), zpts);
       [xsp,ysp,zsp] = ndgrid(xsp,ysp,zsp);
       val = [xsp(:),ysp(:),zsp(:), maxh+0*xsp(:)];
       tmpname = write_tmp_mszfile( val );
       opt.options.meshsizefilename = tmpname;
   otherwise
       eval(sprintf('opt.%s = val;',varargin{idx}));
   end
end

function fname = write_tmp_mszfile( mszpoints )
   % From Documentation: Syntax is
   % np
   % x1 y1 z1 h1
   % x2 y2 z2 h2
   n_pts_elecs=  size(mszpoints,1);
   fname = [tempname,'.msz'];
   fname = strrep(fname,'\','/'); % needs unix-style path on windows
   fid=fopen(fname,'w');
   fprintf(fid,'%d\n',n_pts_elecs);
   for i = 1:size(mszpoints,1)
      fprintf(fid,'%10f %10f %10f %10f\n', mszpoints(i,:));
   end
   fprintf(fid,'0\n'); % lines
   fclose(fid); % ptsfn

% copy all fields from usr to opt
% check that the fields exist in opt
function opt = copy_field(opt,usr)
optnms = fieldnames(opt);
usrnms = fieldnames(usr);
for i = 1:length(usrnms)
   % check if field exist
   if ~ismember(usrnms{i},optnms)
     error('Unsupported option %s',usrnms{i});
   end
   if isstruct(usr.(usrnms{i})) % recurse
      opt.(usrnms{i}) = copy_field(opt.(usrnms{i}), usr.(usrnms{i}) );
   else
      opt.(usrnms{i}) = usr.(usrnms{i});
   end
end


% write the ng.opt file
function write_ng_file(opt)
fid = fopen('ng.opt','w');
flds = fieldnames(opt);
write_field(fid,opt,[]);
fclose(fid);


% recurses over fields and prints to file
function write_field(fid,s,str)
if isstruct(s)
   if ~isempty(str)
      str = [str '.'];
   end
   flds = fieldnames(s);
   for i = 1:length(flds)
      write_field(fid,s.(flds{i}),[str flds{i}]);
   end
elseif ischar(s)
   fprintf(fid,'%s  %s\n', str, s);
elseif round(s) == s
   fprintf(fid,'%s  %d\n', str, s);
else
   fprintf(fid,'%s  %f\n', str, s);
end


function opt = default_opt(varargin)
if nargin == 0
    str = 'moderate';
else
    str = varargin{1};
end
switch str
    case 'very_coarse'
        v = [1, 1, 0.3, 0.7, 0.25, 0.8, 0.2, 0.5, 0.25, 1];
    case 'coarse'
        v = [2, 1.5, 0.5, 0.5, 1, 0.35, 1, 0.5, 1.5];
    case 'moderate'
        v = [3, 2, 1, 0.3, 1, 1.5, 0.5, 2, 1, 2];
    case 'fine'
        v = [4, 3, 2, 0.2, 1.5, 2, 1.5, 3.5, 1.5, 3];
    case 'very_fine'
        v = [5, 5, 3, 0.1, 3, 5, 3, 5, 3, 5];
end

opt.dirname = '.';
opt.loadgeomtypevar = '"All Geometry types (*.stl,*.stlb,*.step,*.stp,...)"';
opt.exportfiletype = '"Neutral Format"';
opt.meshoptions.fineness = v(1);
opt.meshoptions.firststep = 'ag';
opt.meshoptions.laststep = 'ov';
opt.options.localh = 1;
opt.options.delaunay = 1;
opt.options.checkoverlap = 1;
opt.options.checkchartboundary = 1;
opt.options.startinsurface = 0;
opt.options.blockfill = 1;
opt.options.debugmode = 0;
opt.options.dooptimize = 1;
opt.options.parthread = 1;
opt.options.elsizeweight = 0.2;
opt.options.secondorder = 0;
opt.options.elementorder = 1;
opt.options.quad = 0;
opt.options.inverttets = 0;
opt.options.inverttrigs = 0;
opt.options.autozrefine = 0;
opt.options.meshsize = 1000;
opt.options.minmeshsize = 0;
opt.options.curvaturesafety = v(2);
opt.options.segmentsperedge = v(3);
opt.options.meshsizefilename = '';
opt.options.badellimit = 175;
opt.options.optsteps2d = 3;
opt.options.optsteps3d = 5;
opt.options.opterrpow = 2;
opt.options.grading = v(4);
opt.options.printmsg = 2;
opt.geooptions.drawcsg = 1;
opt.geooptions.detail = 0.001;
opt.geooptions.accuracy = 1e-6;
opt.geooptions.facets = 20;
opt.geooptions.minx = -1000;
opt.geooptions.miny = -1000;
opt.geooptions.minz = -1000;
opt.geooptions.maxx = 1000;
opt.geooptions.maxy = 1000;
opt.geooptions.maxz = 1000;
opt.viewoptions.specpointvlen = 0.3;
opt.viewoptions.light.amb = 0.3;
opt.viewoptions.light.diff = 0.7;
opt.viewoptions.light.spec = 1;
opt.viewoptions.light.locviewer = 0;
opt.viewoptions.mat.shininess = 50;
opt.viewoptions.mat.transp = 0.3;
opt.viewoptions.colormeshsize = 0;
opt.viewoptions.whitebackground = 1;
opt.viewoptions.drawcolorbar = 1;
opt.viewoptions.drawcoordinatecross = 1;
opt.viewoptions.drawnetgenlogo = 1;
opt.viewoptions.stereo = 0;
opt.viewoptions.drawfilledtrigs = 1;
opt.viewoptions.drawedges = 0;
opt.viewoptions.drawbadels = 0;
opt.viewoptions.centerpoint = 0;
opt.viewoptions.drawelement = 0;
opt.viewoptions.drawoutline = 1;
opt.viewoptions.drawtets = 0;
opt.viewoptions.drawprisms = 0;
opt.viewoptions.drawpyramids = 0;
opt.viewoptions.drawhexes = 0;
opt.viewoptions.drawidentified = 0;
opt.viewoptions.drawpointnumbers = 0;
opt.viewoptions.drawededges = 1;
opt.viewoptions.drawedpoints = 1;
opt.viewoptions.drawedpointnrs = 0;
opt.viewoptions.drawedtangents = 0;
opt.viewoptions.shrink = 1;
opt.stloptions.showtrias = 0;
opt.stloptions.showfilledtrias = 1;
opt.stloptions.showedges = 1;
opt.stloptions.showmarktrias = 0;
opt.stloptions.showactivechart = 0;
opt.stloptions.yangle = 10;
opt.stloptions.contyangle = 20;
opt.stloptions.edgecornerangle = 0;
opt.stloptions.chartangle = 0;
opt.stloptions.outerchartangle = 120;
opt.stloptions.usesearchtree = 0;
opt.stloptions.chartnumber = 1;
opt.stloptions.charttrignumber = 1;
opt.stloptions.chartnumberoffset = 0;
opt.stloptions.atlasminh = 0.1;
opt.stloptions.resthsurfcurvfac = v(5);
opt.stloptions.resthsurfcurvenable = 0;
opt.stloptions.resthatlasfac = 2;
opt.stloptions.resthatlasenable = 1;
opt.stloptions.resthchartdistfac = v(6);
opt.stloptions.resthchartdistenable = 0;
opt.stloptions.resthlinelengthfac = v(7);
opt.stloptions.resthlinelengthenable = 1;
opt.stloptions.resthcloseedgefac = v(8);
opt.stloptions.resthcloseedgeenable = 1;
opt.stloptions.resthedgeanglefac = v(9);
opt.stloptions.resthedgeangleenable = 0;
opt.stloptions.resthsurfmeshcurvfac = v(10);
opt.stloptions.resthsurfmeshcurvenable = 0;
opt.stloptions.recalchopt = 1;
opt.visoptions.subdivisions = 1;

function do_unit_test
opt.meshoptions.fineness = 6;
ng_write_opt(opt);
fid = fopen('ng.opt','r'); str= fread(fid,[1,inf],'uint8=>char'); fclose(fid);
unit_test_cmp('fineness=6',isempty( ...
    strfind(str, 'meshoptions.fineness  6')), 0);

ng_write_opt('meshoptions.fineness',4,'meshoptions.firststep','aaa');
fid = fopen('ng.opt','r'); str= fread(fid,[1,inf],'uint8=>char'); fclose(fid);
unit_test_cmp('firststep=aaa',isempty( ...
    strfind(str, 'meshoptions.firststep  aaa')), 0);

% the last one will not be written since output was requested
opt = ng_write_opt('meshoptions.fineness',4,'meshoptions.firststep','bbb');
fid = fopen('ng.opt','r'); str= fread(fid,[1,inf],'uint8=>char'); fclose(fid);
unit_test_cmp('firststep=bbb',isempty( ...
    strfind(str, 'meshoptions.firststep  bbb')), 1);
