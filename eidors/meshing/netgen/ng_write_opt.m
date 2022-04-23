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
%   or
%  ng_write_opt('MSZSPHERE', [xctr, yctr, zctr, radius, maxh])
%  ng_write_opt('MSZCYLINDER', [x1, y1, z1, x2, y2, z2, radius, maxh])
%     where (x1,y1,z1) and (x2,y2,z2) are the limits on the axis
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
if nargs >= 1 && ischar(varargin{1}) && ~any(varargin{1} == '.')  ...
              && isempty(strfind(varargin{1},'MSZ'))
   str = varargin(1);
   varargin(1) = [];
   nargs = nargs - 1;
end

% get default options
opt = default_opt(str{:});

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
       val = msz_brick(val);
       tmpname = write_tmp_mszfile( val );
       opt.options.meshsizefilename = tmpname;
   case 'MSZSPHERE'
       val = msz_sphere(val);
       tmpname = write_tmp_mszfile( val );
       opt.options.meshsizefilename = tmpname;
   case 'MSZCYLINDER'
       val = msz_cylinder(val);
       tmpname = write_tmp_mszfile( val );
       opt.options.meshsizefilename = tmpname;
   otherwise
       eval(sprintf('opt.%s = val;',varargin{idx}));
   end
end

function val = msz_brick(val)
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

function val = msz_sphere(val)
%  ng_write_opt('MSZSPHERE', [xctr, yctr, zctr, radius, maxh])
    maxh = val(5);
    radius = val(4);
    npts= floor(2*radius/maxh)+1;
    xsp= linspace(val(1) - radius, val(1) + radius, npts);
    ysp= linspace(val(2) - radius, val(2) + radius, npts);
    zsp= linspace(val(3) - radius, val(3) + radius, npts);
    [xsp,ysp,zsp] = ndgrid(xsp,ysp,zsp);
    s_idx = ((xsp-val(1)).^2 + ...
             (ysp-val(2)).^2 + ...
             (zsp-val(3)).^2) < radius^2;
    s_idx = s_idx(:);
    val = [xsp(s_idx),ysp(s_idx),zsp(s_idx), maxh+0*xsp(s_idx)];

% space points around zero with spacing maxh
function pts = space_around_zero(lim,maxh);
    pts = 0:maxh:lim;
    pts = [-fliplr(pts(2:end)),pts];

function val = msz_cylinder(val)
    % [x1, y1, z1, x2, y2, z2, radius, maxh])
    len2= norm(val(1:3) - val(4:6))/2;
    ctr =     (val(1:3) + val(4:6))/2;
    maxh = val(8);
    radius = val(7);
    xsp= space_around_zero(radius, maxh);
    ysp= space_around_zero(radius, maxh);
    zsp= space_around_zero(len2, maxh);
    [xsp,ysp,zsp] = ndgrid(xsp,ysp,zsp);
    s_idx = (xsp.^2 + ysp.^2) < radius^2;
    xsp = xsp(s_idx(:));
    ysp = ysp(s_idx(:));
    zsp = zsp(s_idx(:));
% Rotate: 
% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/897677#897677
    ab = [0;0;2*len2] + [val(4:6)-val(1:3)]';
    R=2*ab*ab'/(ab'*ab) - eye(3);
    val = [xsp,ysp,zsp]*R + ctr;
    val(:,4) = maxh;

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
%  for i = 1:size(mszpoints,1)
%     fprintf(fid,'%10f %10f %10f %10f\n', mszpoints(i,:));
%  end
%  vectorize to speed up
   fprintf(fid,'%10f %10f %10f %10f\n', mszpoints');
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
   if isempty(str)
      str = '';
   else
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
        % fineness 1
        v = [6, 1.0, 0.3, 0.7, 0.25, 0.8, 0.20, 0.5, 0.25, 1.0];
    case 'coarse'
        % fineness 2
        v = [6, 1.5, 0.5, 0.5, 0.50, 1.0, 0.35, 1.0, 0.50, 1.5];
    case 'moderate'
        % fineness 3
        v = [6, 2.0, 1.0, 0.3, 1.00, 1.5, 0.50, 2.0, 1.00, 2.0];
    case 'fine'
        % fineness 4
        v = [6, 3.0, 2.0, 0.2, 1.50, 2.0, 1.50, 3.5, 1.50, 3.0];
    case 'very_fine'
        % fineness 5
        v = [6, 5.0, 3.0, 0.1, 3.00, 5.0, 3.00, 5.0, 3.00, 5.0];
    otherwise % use moderate
        % this is called if some other ng_write_opt is used, like MSZBRICK
        v = [6, 2.0, 1.0, 0.3, 1.00, 1.5, 0.50, 2.0, 1.00, 2.0];
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
opt.stloptions.yangle = 30;
opt.stloptions.contyangle = 20;
opt.stloptions.edgecornerangle = 60;
opt.stloptions.chartangle = 15;
opt.stloptions.outerchartangle = 70;
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
%  test_main_options()
   unit_test_MSZ()

function test_main_options()
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

function unit_test_MSZ
   shape_str =  ...
     'solid mainobj = orthobrick(-9,-9,-9;9,9,9);';
   i=0; while(true); i=i+1; switch i
       case 1; n_exp = 8;
           ng_write_opt('very_coarse');
       case 2; n_exp = 22;
           ng_write_opt('fine');
       case 3; n_exp = 59;
           ng_write_opt('very_fine');
       case 4; n_exp = 8;
           ng_write_opt('moderate');
       case 5; n_exp = 746;
           ng_write_opt('MSZSPHERE',[5,5,5,4,1]);
       case 6; n_exp = 106;
           ng_write_opt('MSZBRICK',[5,5,5,9,9,9,1]);
       case 7; n_exp = 544;
           ng_write_opt('MSZCYLINDER',[5,5,5,8,8,8,5,1]);
       case 8; n_exp = 23;
           ng_write_opt('MSZCYLINDER',[5,5,5,8,8,8,5,5]);
         
       otherwise; break
       end
       eidors_cache clear
       fmdl = ng_mk_gen_models(shape_str,[],[],'');
       unit_test_cmp('opt:', num_nodes(fmdl), n_exp)
%      show_fem(fmdl); disp(num_nodes(fmdl))
   end
