function opt = ng_write_opt(varargin)
%NG_WRITE_OPT Write an ng.opt file in current directory
%  NG_WRITE_OPT, without inputs, creates a default ng.opt
%
%  NG_WRITE_OPT(OPTSTRUCT) creates uses options as specified in OPTSTRUCT
%
%  NG_WRITE_OPT(PATH) copies the file specified in PATH as ng.opt to the
%  current directory.
%
%  NG_WRITE_OPT('sec.option',val,...) offers an alternative interface to
%  modify specific options
%
%  OPTSTRUCT = NG_WRITE_OPT(...) returns a struct with the options.
%  No file is written.
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
%  Caveat: Currently it seems that Netgen on Windows ignores the ng.opt
%  file
%
%  See also CALL_NETGEN

% (C) 2012 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

%TODO:
% 1. Check which options require meshoptions.fineness = 6 and enforce it

% if input is 'UNIT_TEST', run tests
if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'UNIT_TEST') 
   do_unit_test; return; end

if nargin == 1 && ischar(varargin{1}) % a path
   copyfile(varargin{1},'ng.opt');
   return
end
% get default options
opt = default_opt;
if nargin == 0
    usr = struct;
end
% modify as per user input
if nargin == 1 && isstruct(varargin{1})
   usr = varargin{1};
end
if nargin > 1 % string value pairs
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
   eval(sprintf('opt.%s = val;',varargin{idx}));
end


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


function opt = default_opt
opt.dirname = '.';
opt.loadgeomtypevar = '"All Geometry types (*.stl,*.stlb,*.step,*.stp,...)"';
opt.exportfiletype = '"Neutral Format"';
opt.meshoptions.fineness = 3;
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
opt.options.curvaturesafety = 0.2;
opt.options.segmentsperedge = 1;
opt.options.meshsizefilename = '';
opt.options.badellimit = 175;
opt.options.optsteps2d = 3;
opt.options.optsteps3d = 5;
opt.options.opterrpow = 2;
opt.options.grading = 0.3;
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
opt.stloptions.resthsurfcurvfac = 1;
opt.stloptions.resthsurfcurvenable = 0;
opt.stloptions.resthatlasfac = 2;
opt.stloptions.resthatlasenable = 1;
opt.stloptions.resthchartdistfac = 1.5;
opt.stloptions.resthchartdistenable = 0;
opt.stloptions.resthlinelengthfac = 0.5;
opt.stloptions.resthlinelengthenable = 1;
opt.stloptions.resthcloseedgefac = 2;
opt.stloptions.resthcloseedgeenable = 1;
opt.stloptions.resthedgeanglefac = 1;
opt.stloptions.resthedgeangleenable = 0;
opt.stloptions.resthsurfmeshcurvfac = 2;
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
