function imdl = mk_geophysics_model(str, ne, opt);
% imdl = mk_geophysics_model(str, ne, [option])
%
% ne  - number of electrodes, 5 metre spacing (+5,+10,...)
%       and 0.1 metre diameter
%         OR
%       a list of electrode locations in the x-dimension or a 2- or
%       3-dimensional array, one electrode per row, missing columns
%       will be set to zero
%       ne = 16
%       ne = [4 6 10 20] % 1d: (x)
%       ne = [0 1; 2 1.5; 3 1.2; 7 2.5] % 2d: (x,y)
%       ne = [0 0.1 1; 2 -0.1 1.5; 3 -0.15 1.2; 7 0 2.5] % 3d: (x,y,z)
% str - model, x = see hmax_rec
%       h2x -   2D half-space, linear CEM array (2d fwd)
%       h2p5x - 2.5D half-space, linear CEM array (2d fwd + Fourier y-dimension)
%       h3x   - 3D half-space, linear CEM array (3d fwd)
%       h22x  - 2D half-space, linear CEM array (2d fwd, 2d rec)
%       h32x  - 2D half-space, linear CEM array (3d fwd, 2d rec)
%       h33x  - 2D half-space, linear CEM array (3d fwd, 3d rec)
%       H2x -   2D half-space, linear PEM array (2d fwd)
%       H2p5x - 2.5D half-space, linear PEM array (2d fwd + Fourier y-dimension)
%       H3x   - 3D half-space, linear PEM array (3d fwd)
%       H22x  - 2D half-space, linear PEM array (2d fwd, 2d rec)
%       H32x  - 2D half-space, linear PEM array (3d fwd, 2d rec)
%       H33x  - 2D half-space, linear PEM array (3d fwd, 3d rec)
% opt - override default configuration options (optional cell array)
%       'hmax_fwd' - fine reconstruction mesh density, given an
%                    array width xw, and electrode spacing es
%                    Note that for ne=16, 'A' and 'a' are equivalent.
%                ['a' : hmax_fwd=xw/1;    ['A' : hmax_fwd=es*16;
%                 'b' : hmax_fwd=xw/2;     'B' : hmax_fwd=es*8;
%                 'c' : hmax_fwd=xw/4;     'C' : hmax_fwd=es*4;
%                 ...                       ...
%                 'z' : hmax_fwd=xw/2^25]  'Z' : hmax_fwd=es*2^-21]
%       'hmax_fwd_inner'
%                  - reconstruction model mesh density for the inner region
%                    [default: 1/2 of the outer region density hmax_fwd]
%       'hmax_rec' - reconstruction model mesh density [hmax_fwd*2]
%       'elec_width' - electrode width [0.1 m]
%                    width = 0 requests a PEM, rather than CEM
%       'z_contact' - electrode contact impedance [0.01 \Ohm.m]
%       'elec_spacing' - distance between electrode centers [5 m]
%       'extend_x' - extra mesh in the principle axis of the
%                    electrode array, multiple of array width [1]
%       'extend_y' - extra mesh in the minor axis of the
%                    electrode array, multiple of array width
%                    (3D models only) [1]
%       'extend_z' - extra depth of model, multiple of array width [1]
%       'extend_inner_x'
%                  - inner (denser) mesh, multiple of array width [3/5]
%       'extend_inner_y'
%                  - inner (denser) mesh, multiple of array width [3/5]
%                    (3D models only)
%       'extend_inner_z'
%                  - inner (denser) mesh, multiple of array width [2/5]
%       'skip_c2f' - skip building the rec_model to fwd_model mapping [0]
%       'threshold' - threshold for electrode placement errors [1e-12]
%                    if error exceeds threshold, mash nodes by cubic interp
%                    maintaining the side and lower boundaries of the
%                    model... the nodes will be mashed until the electrodes
%                    conform to the requested electrode positions,
%                    correcting for Netgen inaccuracies and lack of
%                    vertical profile (Inf = disable node mashing)
%       'build_stim' - use stim_pattern_geophys to build a standard geophysics
%                      stim/meas sequence and add it to the model, based on the
%                      number of electrodes, and assuming a co-linear array;
%                      set to 'none' to skip this step
%                      [default: 'Wenner']
%       'extra_ng_code' - pass extra netgen code, any 'tlo' are subtracted from
%                      the inner region, units are scaled, so left to right-most
%                      electrodes are -1 to +1
%                      [default: {}]
%
% The linear electrode array runs in the +X direction at Z=0. For
% the 3D model, the Y-axis is perpendicular to the electrode array.
%
% (C) 2015--2017 A. Boyle
% License: GPL version 2 or version 3

% model: 64 electrode, 2d half-space
% Once upon a time, this code started out from the following tutorial.
% model from http://eidors3d.sourceforge.net/tutorial/other_models/square_mesh.shtml

if ischar(str) && strcmp(str,'UNIT_TEST'); do_unit_test; return; end
copt.fstr = 'mk_geophysics_model';
if nargin < 3
   opt = {};
end
SALT='z$Id$'; % stick a key in the model 'save' file, so we can expire them when the model definitions age
imdl = eidors_cache(@mk_model,{str, ne, opt, SALT}, copt);

function imdl=mk_model(str,ne,opt,SALT);
if str(1) ~= 'h' && str(1) ~= 'H'
   error([str ': I only know how to build linear half-space models: h***']);
end

MDL_2p5D_CONFIG = 0;
switch str(2:end-1)
   case {'2', '3'} % simple meshes
      FMDL_DIM = str(2) - '0';
      CMDL_DIM = 0; % no cmdl
   case '2p5' % 2.5D Fourier transformed
      FMDL_DIM = 2;
      CMDL_DIM = 0; % no cmdl
      MDL_2p5D_CONFIG = 1;
   case {'22', '33', '32'} % dual meshes
      FMDL_DIM = str(2) - '0';
      CMDL_DIM = str(3) - '0';
   otherwise
      error([str ': unrecognized model type']);
end
assert(CMDL_DIM ~= 3, '3d rec_model not yet tested');

skip_c2f = 0;
if str(1) == 'h'
   elec_width = 0.1;
else % str(1) == 'H'
   elec_width = 0;
end
z_contact= 0.01;
nodes_per_elec= 3; %floor(elec_width/hmax_rec*10);
elec_spacing= 5.0;
threshold = 1e-12;
save_model_to_disk = (length(ne) == 1) && (length(opt) == 0);

extend_x = 1;
extend_y = 1;
extend_z = 1;
extend_inner_x = 3/5;
extend_inner_y = 3/5;
extend_inner_z = 2/5;
build_stim = 'Wenner';
extra_ng_code = '';
if length(opt) > 0 % allow overriding the default values
   assert(round(length(opt)/2)*2 == length(opt),'option missing value?');
   expect = {'hmax_rec','hmax_fwd', 'hmax_fwd_inner', ...
             'elec_width','z_contact','elec_spacing',...
             'extend_x', 'extend_y', 'extend_z', ...
             'extend_inner_x', 'extend_inner_y', 'extend_inner_z', ...
             'skip_c2f', 'threshold', 'build_stim', 'extra_ng_code'};
   opts = struct(opt{:})
   for i = fieldnames(opts)'
      assert(any(strcmp(i,expect)), ['unexpected option: ',i{:}]);
      eval([i{:} ' = opts.(i{:});']);
   end
   if (str(1) == 'H') && isfield(opts, 'elec_width')
      error('requested "H" PEM model but configured "elec_width" option');
   end
end
if length(ne) == 1 % ne: number of electrodes
   xw=(ne-1)*elec_spacing; % array width
   %xs=-(ne-1)*elec_spacing/2; % array centered
   xs=+5; % array at left-most at +5
   xyz = xs+([1:ne]'-1)*elec_spacing;
else
   xyz = ne; % must be a set of coordinates for the electrodes...
end
if size(xyz,1) == 1
   xyz = xyz'; % flip to column
end
xyz = [xyz zeros(size(xyz,1),3-size(xyz,2))]; % [x 0 0] or [ x y 0 ] or [ x y z ]
ne=size(xyz,1);
[R, X] = rot_line_to_xaxis(xyz);
% % rescale, centre electrodes so NetGen can be happy
xyzc = (xyz - X)*R; % centre and scale electrodes: -1 to +1 y-axis
% xyzc = xyz; % centre and scale electrodes: -1 to +1 y-axis
xw=max(xyzc(:,1))-min(xyzc(:,1));
xs=min(xyzc(:,1));
elec_spacing = min(min(pdist(xyzc) + diag(inf*(1:size(xyzc,1))))); % min spacing btw elec

if ~exist('hmax_fwd','var')
   if str(end)-'a' >= 0
      hmax_fwd = xw*2^-(str(end)-'a');
   else
      hmax_fwd = elec_spacing*2^-(str(end)-'A'-4);
   end
end
if ~exist('hmax_rec','var') % allow hmax_rec to depend on configured hmax_fwd
   hmax_rec=hmax_fwd*2.01; % avoid parametrization aliasing
end
if ~exist('hmax_fwd_inner','var') % allow hmax_rec to depend on configured hmax_fwd
   hmax_fwd_inner=hmax_fwd/2.0;
end

if save_model_to_disk
   isOctave = exist('OCTAVE_VERSION');
   if isOctave
      octavestr='-octave-';
   else
      octavestr='-';
   end
   % NOTE models are stored in the directory specified by eidors_cache('cache_path')
   filename=fullfile(eidors_cache('cache_path'),sprintf('mk_geophysics_model%simdl-%s-%03del.mat',octavestr,str,ne));
   if exist(filename, 'file') == 2
      tmp = whos('-file',filename,'SALT');
      if length(tmp) > 0
         tmp = load(filename,'SALT');
         fSALT = tmp.SALT;
      else
         fSALT = 'deadbeef';
      end
      if strcmp(fSALT, SALT)
         tmp=load(filename,'imdl');
         imdl = tmp.imdl;
         eidors_msg(sprintf('%s: %s, %d electrode model loaded from file',filename,str,ne));
         return
      end % hmm, the SALT doesn't match so we go back to generating a new model from scratch
   end
end

assert(extend_x>0,'extend_x must be > 0');
assert(extend_y>0,'extend_y must be > 0');
assert(extend_z>-1,'extend_z must be > -1');

% Calculate cmdl, fmdl and inner fmdl (fmdlin) min/max coordinates.
% After rescaling (normalizing), electrode major axis is along the x-axis, and
% the electrode array will be 2 units long (-1,+1).
fmdlbox =   [-(1+2*extend_x)        +(1+2*extend_x);
             -(1+2*extend_y)        +(1+2*extend_y);
             -(2+2*extend_z)        +3];
fmdlinbox = [-(1+2*extend_inner_x)  +(1+2*extend_inner_x);
             -(1+2*extend_inner_y)  +(1+2*extend_inner_y);
             -(2+2*extend_inner_z)  +2];
cmdlbox = fmdlbox;

assert(all(fmdlbox(:,1) <= fmdlinbox(:,1)), 'oops, inner mesh must be smaller than outer mesh');
assert(all(fmdlbox(:,2) >= fmdlinbox(:,2)), 'oops, inner mesh must be smaller than outer mesh');
assert(all(fmdlbox(:,1) <= cmdlbox(:,1)), 'oops, reconstruction mesh must be smaller than forward mesh');
assert(all(fmdlbox(:,2) >= cmdlbox(:,2)), 'oops, reconstruction mesh must be smaller than forward mesh');

% build shape string for NetGen
if elec_width ~= 0
   elec_shape = [elec_width*norm(R)/2, 0, elec_width*norm(R)/2/(nodes_per_elec-1)]; % CEM, circular, maxh
   elec_pos   = [ xyzc(:,1:FMDL_DIM), repmat([zeros(1,3-FMDL_DIM+2) 1],ne,1) ]; % p(x,y,z=0), n(0,0,1)
   cem2pem = 0;
else
   if FMDL_DIM == 3
      elec_width = elec_spacing/2;
      elec_shape = [elec_width, elec_width, hmax_fwd_inner]; % PEM, sz, maxh
      elec_pos   = [ xyzc, repmat([0 0 1],ne,1) ]; % p(x,y,z=0), n(0,0,1)
      elec_pos(:,1) = elec_pos(:,1) + elec_width/2;
      elec_pos(:,2) = elec_pos(:,2) + elec_width/2;
      cem2pem = 1;
   else % FMDL_DIM == 2
      elec_width = elec_spacing/2;
      elec_shape = [elec_width, elec_width, hmax_fwd_inner]; % rectangular CEM->PEM, maxh
      elec_pos   = [ xyzc(:,1:2), repmat([0 0 0 1],ne,1) ]; % p(x,y,z=0), n(0,0,1)
      elec_pos(:,1) = elec_pos(:,1) + elec_width/2;
      cem2pem = 1;
   end
end
skip_tlo = '';
for idx = strfind(extra_ng_code, 'tlo ')
   sc = strfind(extra_ng_code, ';'); % semicolons
   sc = min(sc + (sc<idx)*1e10);
   tlo = extra_ng_code(idx+4:sc-1);
   skip_tlo = [skip_tlo ' and (not ' tlo ')'];
end
elec_pos(:,FMDL_DIM) = 0; % flatten electrode positions onto the "ps" plane
if FMDL_DIM == 3
   shape_str = [...
                sprintf('solid ps = plane(%s);\n', a2s([0 0 0; 0 0 1])), ...
                sprintf('solid bi = orthobrick(%s);\n', a2s(fmdlinbox')), ...
                sprintf('solid bo = orthobrick(%s);\n', a2s(fmdlbox')), ...
                extra_ng_code, ...
                sprintf('solid ri = bi and ps%s -maxh=%f;\n', skip_tlo, hmax_fwd_inner), ...
                sprintf('solid ro = bo and ps and (not bi) -maxh=%f;\n', hmax_fwd), ...
                sprintf('solid mainobj = ri;\n')];
   % Note that ri must be the 'mainobj' so that it can intersect with the electrodes
   % additional top level objects for netgen
else % netgen 2d model
   % to create the 2D slice we need to give NetGen something to work with
   %Need some width to let netgen work, but not too much so
   % that it meshes the entire region
   sw = range(fmdlinbox(1,:)) / 5; % initial estimate
   sw = min(sw,2*hmax_fwd); % coarse model slice width
   ri2d = fmdlinbox'; ri2d(:,2) = [-sw 0 ];
   ro2d = fmdlbox'; ro2d(:,2) = [-sw 0 ];
   shape_str = [...
                sprintf('solid ps = plane(%s);\n', a2s([0 0 0; 0 0 1])), ...
                sprintf('solid bi = orthobrick(%s);\n', a2s(ri2d)), ...
                sprintf('solid bo = orthobrick(%s);\n', a2s(ro2d)), ...
                extra_ng_code, ...
                sprintf('solid ri = bi and ps%s -maxh=%f;\n', skip_tlo, hmax_fwd_inner), ...
                sprintf('solid ro = bo and ps and (not bi) -maxh=%f;\n', hmax_fwd), ...
                sprintf('solid mainobj = ri;\n')];
end
% fprintf('SHAPE_STR: %s', shape_str); elec_pos
elec_obj = 'ps';
[fmdl, mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj, 'tlo ro');
if FMDL_DIM == 2 % 2D
   % now convert the roughly 2D slice into a true 2D plane
   [fmdl, mat_idx] = copy_mdl2d_from3d(fmdl, mat_idx, 'y');
else % 3D
   if CMDL_DIM ~= 0
      % c2f
      cmdl.mk_coarse_fine_mapping.f2c_offset  = [0 0 0];
      cmdl.mk_coarse_fine_mapping.f2c_project = [1 0 0; 0 0 1; 0 1 0];
      % duplicate parameters since mk_analytic/approx_c2f have different names...
      cmdl.mk_analytic_c2f.f2c_offset  = cmdl.mk_coarse_fine_mapping.f2c_offset;
      cmdl.mk_analytic_c2f.f2c_project = cmdl.mk_coarse_fine_mapping.f2c_project;
   end
end

if cem2pem
   fmdl = convert_cem2pem(fmdl, xyzc);
end

% 2d cmdl
xllim = fmdlbox(1,1);
xrlim = fmdlbox(1,2);
zdepth = fmdlbox(3,1);
% AA - 2022-05 Updated from hmax_rec/2 
xr=max(floor((xrlim-xllim)/hmax_rec),1)*2+1; % odd number
yr=max(floor(-zdepth/hmax_rec),1)*2+1; % odd number
[x,y] = meshgrid( linspace(xllim,xrlim,xr), linspace(zdepth,0,yr));
vtx= [x(:),y(:)];
if CMDL_DIM ~= 0
   cmdl= mk_fmdl_from_nodes( vtx,{vtx(1,:)}, z_contact, 'sq_m2');
end

% stick electrode nodes into cmdl so that show_fem will plot them
for i=1:ne
   n=fmdl.electrode(i).nodes;
   nn=length(n);
   nx=fmdl.nodes(n,:);

   fmdl.electrode(i).z_contact = z_contact;
   if CMDL_DIM ~= 0
      nnc = length(cmdl.nodes);
      cmdl.nodes = [cmdl.nodes; nx(:,[1 FMDL_DIM])];
      cmdl.electrode(i).nodes = (nnc+1):(nnc+nn);
      cmdl.electrode(i).z_contact = z_contact;
   end
end

% fix electrode locations if necessary
elec_err = sqrt(sum(mdl_elec_err(fmdl, xyzc).^2,2));
if max(elec_err) > threshold % put electrodes in the right place
   [fmdl, cf] = correct_electrode_positions(fmdl, xyzc);

   if CMDL_DIM ~= 0
      [cmdl, cc] = correct_electrode_positions(cmdl, xyzc);
   end
end
% save functions for later use
fmdl.mv_elec = @shift_electrode_positions;

% reverse the centre and scaling
nn = size(fmdl.nodes,1);
Xn = repmat(X(1,:), nn, 1);
fmdl.nodes = ([fmdl.nodes zeros(nn,3-FMDL_DIM)]/ R) + Xn;
fmdl.nodes = fmdl.nodes(:,1:FMDL_DIM);

if CMDL_DIM ~= 0
   nn = size(cmdl.nodes,1);
   Xn = repmat(X(1,:), nn, 1);
   if CMDL_DIM ~= FMDL_DIM % 2D
      cmdl.nodes = ([cmdl.nodes(:,1) zeros(nn,1) cmdl.nodes(:,2:end)]/ R) + Xn;
      cmdl.nodes = cmdl.nodes(:,[1 3]);
%     cmdl.nodes = cmdl.nodes(:,[1 3])/R([1 3],[1 3]); % ONLY DIVIDE ONCE
   else % CMDL_DIM == FMDL_DIM
      cmdl.nodes = ([cmdl.nodes zeros(nn,3-CMDL_DIM)]/ R) + Xn;
      cmdl.nodes = cmdl.nodes(:,1:CMDL_DIM);
   end
end

% check that all electrodes were found
for i = 1:length(fmdl.electrode)
   nn = fmdl.electrode(i).nodes;
   assert(length(nn(:)) > 0, sprintf('electrode#%d: failed to find nodes',i));
end

% show_fem(fmdl); xlabel('x'); ylabel('y'); zlabel('z');

% TODO ... Actually we want a node in the middle: the boundary is bad too...
[~, gn] = min(fmdl.nodes(:,end));
fmdl.gnd_node = gn; % make sure the ground node is away from surface electrodes
if CMDL_DIM ~= 0
   [~, gn] = min(cmdl.nodes(:,end));
   cmdl.gnd_node = gn; % make sure the ground node is away from surface electrodes
end


imdl= mk_common_model('a2d0c',4); % 2d model
imdl.fwd_model = fmdl;
imdl.name = ['EIDORS mk_geophysics_model ' str];
imdl.fwd_model.name = ['EIDORS mk_geophysics_model fwd model ' str];
if CMDL_DIM ~= 0
   imdl.rec_model.name = ['EIDORS mk_geophysics_model rec model ' str];
   imdl.rec_model = cmdl;
   % EIDORS "analytic_c2f" gets stuck, do an approximate one
   %eidors_default('set','mk_coarse_fine_mapping','mk_analytic_c2f');
   %eidors_default('set','mk_coarse_fine_mapping','mk_approx_c2f');
   %[c2f,out] = mk_coarse_fine_mapping(fmdl,cmdl);
   if ~skip_c2f
      [c2f,out] = mk_approx_c2f(fmdl,cmdl);
      imdl.fwd_model.coarse2fine = c2f;
      imdl.fwd_model.background = out;
   end
end

if MDL_2p5D_CONFIG
   imdl.fwd_model.jacobian   = @jacobian_adjoint_2p5d_1st_order;
   imdl.fwd_model.solve      = @fwd_solve_2p5d_1st_order;
   imdl.fwd_model.system_mat = @system_mat_2p5d_1st_order;
   imdl.fwd_model.jacobian_adjoint_2p5d_1st_order.k = [0 3];
   imdl.fwd_model.fwd_solve_2p5d_1st_order.k = [0 3];
end

if ~strcmp(build_stim,'none')
   imdl.fwd_model.stimulation = stim_pattern_geophys(length(imdl.fwd_model.electrode), build_stim);
end

if save_model_to_disk
   save(filename,'imdl','SALT');
end

% convert 2x3 array to "x1,y1,z1;x2,y2,z2" string
% convert 1x3 array to "x1,y1,z1" string
% if 1x1 array, then use xy_ctr=x1,y1 and [0 0 1]=x2,y2,z2 (z+)
function s = a2s(a)
if length(a(:)) == 3
   s = sprintf('%f,%f,%f', ...
               a(1), a(2), a(3));
else
   s = sprintf('%f,%f,%f;%f,%f,%f', ...
                a(1,1), a(1,2), a(1,3), ...
               a(2,1), a(2,2), a(2,3));
end

% returns R rotation/scaling and X0 offset
% xyz1 = R * xyz + X; % rotate and scale to +/- 1
function [R,X,var] = rot_line_to_xaxis(xyz)
x = xyz(:,1); y=xyz(:,2); z=xyz(:,3);

% fit line to points
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/32502
p = mean(xyz);
[U,S,V] = svd([x-p(1), y-p(2), z-p(3)]);
if V(end,1) ~= 0
   N=1/V(end,1)*V(:,1);
else
   N=V(:,1);
end
A=p' + dot( xyz(1,  :) - p, N ) * N/norm(N)^2;
B=p' + dot( xyz(end,:) - p, N ) * N/norm(N)^2;

% rotate line to +y-axis
a = N/norm(N);
b = [ 1 0 0 ]'; % +x
% http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
v = cross(a, b);
s = norm(v); c = dot(a, b); % sin, cos
xt = @(v) [   0  -v(3)  v(2); ... % skew symmetric cross-product of v
            v(3)    0  -v(1); ... % diagnoal is the scaling identity matrix
           -v(2)  v(1)    0];
if abs(s) < eps*1e3
   R = eye(3);
else
   R = (eye(3) + xt(v) + xt(v)^2 * (1-c)/s^2);
   R=R'; % for right multiply: xyz*R
end

X = repmat(p, size(xyz,1), 1);
R = R / max(range((xyz-X)*R)) * 2;

DEBUG=0;
if DEBUG
   clf;
   subplot(211);
   plot3(x,y,z,'bo');
   hold on;
   plot3(x(1),y(1),z(1),'go');
   plot3(p(1),p(2),p(3),'ro');
   plot3([A(1),B(1)],[A(2),B(2)],[A(3),B(3)]);
   grid; axis equal; hold off;
   xlabel('x'); ylabel('y'); zlabel('z');
   title('pre-rotation');
   subplot(212);
   xx = (xyz-X)*R;
   plot3(xx(:,1), xx(:,2), xx(:,3),'bo');
   hold on;
   plot3(xx(1,1), xx(1,2), xx(1,3),'go');
   plot3(0,0,0,'ro');
   grid; axis equal; hold off;
   xlabel('x'); ylabel('y'); zlabel('z');
   title('post-rotation');
end

% the (max-min) range of a variable's values
function r = range(a)
r = max(a(:))-min(a(:));

function [mdl2,idx2] = copy_mdl2d_from3d(mdl3,idx3,sel);
% AB: taken from EIDORS function ng_mk_gen_models() subfunction of the same name
% AB: NEW: sel = 'x', 'y' or 'z' -- default was Z, we want X
   if sel == 'x'
      % swap Z and X
      T = [ 0 0 1; 0 1 0; 1 0 0 ];
   elseif sel == 'y'
      % swap Z and Y
      T = [ 1 0 0; 0 0 1; 0 1 0 ];
   elseif sel == 'z'
      T = eye(3);
   else
      error('sel must be "x", "y" or "z"');
   end
   mdl3.nodes = mdl3.nodes * T; % AB: SWAP axes

   % set name
   mdl2 = eidors_obj('fwd_model',sprintf('%s 2D',mdl3.name));

   % set nodes
   [bdy,idx] = find_boundary(mdl3.elems);
   vtx = mdl3.nodes;
   z_vtx = reshape(vtx(bdy,3), size(bdy) );
   z_vtx_thres = max(z_vtx(:))-10*eps*range(z_vtx(:));
   lay0  = find( all(z_vtx >= z_vtx_thres, 2) );
   bdy0  = bdy( lay0, :);

   vtx0  = unique(bdy0(:));
   mdl2.nodes = vtx(vtx0,1:2);

   % set elems
   nmap  = zeros(size(vtx,1),1); nmap(vtx0) = 1:length(vtx0);
   bdy0  = reshape(nmap(bdy0), size(bdy0) ); % renumber to new scheme
   mdl2.elems = bdy0;

   % set boundary
   mdl2.boundary = find_boundary( mdl2.elems);

   % set gnd_node
   mdl2.gnd_node = nmap(mdl3.gnd_node);

   % set material indices
   % TODO: vectorize code
   idx2 = {};
   idx0  = idx( lay0, :);
   for i=1:size(idx3,2)
      idx2{i} = [];
      ii = 1;
      for j=1:size(idx3{i},1)
         idx_tmp = find( idx0==idx3{i}(j) );
         if not(isempty(idx_tmp))
            idx2{i}(ii,1) = idx_tmp(1,1);
            ii = ii + 1;
         end
      end
   end

   % set electrode
   if isfield(mdl3,'electrode')
      mdl2.electrode = mdl3.electrode;
      for i=1:length(mdl2.electrode);
         nn = mdl3.nodes(mdl3.electrode(i).nodes,:);
         enodes = nmap( mdl2.electrode(i).nodes );
         enodes(enodes==0) = []; % Remove 3D layers
         mdl2.electrode(i).nodes = enodes(:)';
      end
   end

   ignore = {'electrode', 'nodes', 'boundary', 'elems', 'gnd_node', 'boundary_numbers', 'mat_idx'};
   for n=fieldnames(mdl3)'
      if ~any(strcmp(n,ignore))
         mdl2.(n{:}) = mdl3.(n{:});
      end
   end

function mdl = convert_cem2pem(mdl, xyzc)
   if ~isfield(mdl, 'electrode')
      return;
   end
   nd = size(mdl.nodes,2); % number of dimensions
   for i=1:length(mdl.electrode)
      en = mdl.electrode(i).nodes;
      nn = mdl.nodes(en,:); % all nodes for this electrode
      if nd == 2
         np = xyzc(i,[1 3]); % true elec location (2d)
      else
         np = xyzc(i,:); % true elec location (3d)
      end
      D = pdist([np; nn]);
      [~,idx] = min(D(1,2:end)); % closest CEM node to ideal PEM location
      mdl.electrode(i).nodes = en(idx);
      % NOTE: used to bump node's location but now we place a corner of the
      % square electrode at precisely the correct location
      %mdl.nodes(en,:) = np; % shift PEM node to true location
   end

function D2 = pdist(X) % row vectors
   if nargin == 2; error('only supports Euclidean distances'); end
   %D2 = bsxfun(@plus, dot(X, X, 1)', dot(Y, Y, 1)) - 2*(X'*Y) % 1d
   D2 = bsxfun(@minus, X, permute(X,[3 2 1]));
   D2 = squeeze(sqrt(sum(D2.^2,2)));

function [mdl, c] = shift_electrode_positions(mdl, dx)
   for i = 1:length(mdl.electrode)
      en = mdl.electrode(i).nodes;
      ex = mdl.nodes(en,:);
      ep(i,:) = (max(ex,[],1) + min(ex,[],1))/2; % mid-point
   end
   [mdl, c] = correct_electrode_positions(mdl, ep + dx);

function [mdl, c] = correct_electrode_positions(mdl, xyzc)
   %eidors_msg('correct_electrode_positions');
   nd = size(mdl.nodes,2);
   c = 0; err = 1;
   while max(err) > eps
      for n = 1:nd
         switch(n)
            case 1
               mdl = mash_nodes(mdl, 'shift_all',     1, 1, xyzc); % X (downslope)
            case 2 % 2d: Y, 3d: Z
               mdl = mash_nodes(mdl, 'shift_surface', 1, nd, xyzc); % Y or Z (vertical)
            case 3 % note: we only do this for 3d
               mdl = mash_nodes(mdl, 'shift_middle',  1, 2, xyzc); % Y (cross-slope)
            otherwise
               error('duh!');
         end
      end
      err = sqrt(sum(mdl_elec_err(mdl, xyzc).^2,2));
      c=c+1;
      if c >= 100
         break;
      end
   end

function   mdl = mash_nodes(mdl, method, idm, dim, elec_true)
   elec_err = mdl_elec_err(mdl, elec_true);
   err = elec_err(:,dim);

   % add borders for electrode positions at
   % the volume boundary and 50% of the edge to electrode distance
   xq = mdl.nodes(:,idm);
   x = [min(xq); ...
        mean([min(xq) min(elec_true(:,idm))]); ...
        elec_true(:,idm);
        mean([max(xq) max(elec_true(:,idm))]); ...
        max(xq)];
   v = [0; 0; err; 0; 0];

   % scale error to match the electrode locations
   interp_method = 'linear';
   if strcmp(method, 'shift_surface')
         interp_method = 'pchip';
   end
   vq = interp1(x, v, xq, interp_method, 'extrap');

   switch method
      case 'shift_all'
         vqs = 1;
      case 'shift_middle'
         yq = mdl.nodes(:,dim);
         yqr = max(yq)-min(yq); % range
         yqm = (max(yq) + min(yq))/2; % middle = (max + min)/2
         vqs = 1-abs(yq-yqm)./(yqr/2); % scale the shift depending on x's distance from midline
      case 'shift_surface' % assumes positive surface
         yq = mdl.nodes(:,dim);
         yqr = max(yq) - min(yq);
         yqm = min(yq); % min
         vqs = abs(yq - yqm)./yqr; % scale by distance from surface
      otherwise
         error(['unrecognized method: ',method]);
   end
   mdl.nodes(:,dim) = mdl.nodes(:,dim) + (vq .* vqs);

% calculate the error in electrode position for the fwd_model
function err = mdl_elec_err(mdl, xyzc)
   if ~isfield(mdl, 'electrode')
      error('electrodes not available on this model, must supply positional error');
   end

   nel=length(mdl.electrode); % number of electrodes
   nd=size(mdl.nodes,2); % number of dimensions

   eu = ones(nel,nd)*NaN; % init
   for i=1:length(mdl.electrode)
      nn = mdl.nodes(mdl.electrode(i).nodes,:); % nodes per electrode
      eu(i,:) = (max(nn,[],1) + min(nn,[],1))./2; % approx centre of each electrode
   end
   err = xyzc(:,1:nd) - eu;

function test_fwd_rec_match(imdl,idx,str)
   fwd = imdl.fwd_model.nodes(:,idx);
   mxf = max(fwd); mnf = min(fwd);
   rec = imdl.rec_model.nodes;
   mxr = max(rec); mnr = min(fwd);
   unit_test_cmp(['fwd-rec match: ',str],[mxf,mnf],[mxr,mnr],10*eps);

function do_unit_test
   ne = 16;
   imdl = mk_geophysics_model('h2p5a', ne);
   imdl.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img = mk_image(imdl.fwd_model, 1);
   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve_halfspace(img);
   vd = fwd_solve(img);
clf; h=plot([vh.meas vd.meas],'o--'); legend('analytic','FEM'); set(gca,'box','off'); set(h,'LineWidth',2);

   imdl1 = mk_geophysics_model('h2a',[1:6]);
   imdl2 = mk_geophysics_model('h2a',[1:6]');
   imdl3 = mk_geophysics_model('h2a',[1:6]'*[1 0]);
   imdl3Hnm2d = mk_geophysics_model('H2a',[1:6],{'threshold',Inf}); % try without mashing nodes, no veritcal geometry... electrodes should be precisely located if the electrodes were correctly placed
   imdl3Hnm3d = mk_geophysics_model('H3a',[1:6],{'threshold',Inf}); % try without mashing nodes, no veritcal geometry... electrodes should be precisely located if the electrodes were correctly placed
   imdl4 = mk_geophysics_model('h2a',[1:6]'*[1 0] + ([1:6]*0+2)'*[0 1]);
   R = @(x) [cosd(x) -sind(x); sind(x) cosd(x)]; % rotation matrix
   X = [0 2];
   imdl5 = mk_geophysics_model('h2a',([1:6]'*[1 0] + ([1:6]*0+1)'*X)*R(-135));
   elec_pos_2d = [1 1; 2 2; 3 1; 4 1.5];
   elec_pos_3d = [1 0 0; 2 0.5 1; 3 -0.5 2.5; 10 0 3];
   imdl2dc = mk_geophysics_model('h2a', elec_pos_2d); % 2D complete electrode model
   imdl3dc = mk_geophysics_model('h3a', elec_pos_3d); % 3D complete electrode model
   imdl2dp = mk_geophysics_model('H2a', elec_pos_2d); % 2D point electrode model
   imdl3dp = mk_geophysics_model('H3a', elec_pos_3d); % 3D point electrode model
   imdl2df = mk_geophysics_model('H2a', elec_pos_2d, {'threshold', Inf}); % 2D point electrode model... without mashing
   imdl3df = mk_geophysics_model('H3a', elec_pos_3d, {'threshold', Inf}); % 3D point electrode model... without mashing

   % dual meshes
   imdlh32 = mk_geophysics_model('h32a', elec_pos_3d);
   imdlh22 = mk_geophysics_model('h22a', elec_pos_2d);
   imdlH32 = mk_geophysics_model('H32a', elec_pos_3d);
   imdlH22 = mk_geophysics_model('H22a', elec_pos_2d);

   % std dual meshes w/ 16 elec
   imdlh32_16 = mk_geophysics_model('h32a', 16);
   test_fwd_rec_match(imdlh32_16,[1,3],'h32a')
   imdlh22_16 = mk_geophysics_model('h22a', 16);
   test_fwd_rec_match(imdlh22_16,[1,2],'h22a')
   imdlH32_16 = mk_geophysics_model('H32a', 16);
   test_fwd_rec_match(imdlH32_16,[1,3],'h22a')
   imdlH22_16 = mk_geophysics_model('H22a', 16);
   test_fwd_rec_match(imdlH22_16,[1,2],'h22a')

if 0
   % Nolwenn's grid
   [yy,xx] = meshgrid(0:3:24, [0 4:2:20*2 20*2+4]); % 9 x 21 grid
   xyz = [xx(:) yy(:) zeros(length(yy(:)),1)];

   % test grid
   % 1. bad electrodes ... ? fixme?
   % 2. mash_nodes breaks for co-located nodes in any particular dimension... need 2d interp?
   [yy,xx]=meshgrid(1:3,1:3); xyz = [xx(:) yy(:) zeros(length(yy(:)),1)]; % 3 x 3 grid
   imdl = mk_geophysics_model('H3a',xyz);
   % TODO add 'g3a' and 'G3a' types?
end

   unit_test_cmp('h2a halfspace vs default TEST', norm(vh.meas - vd.meas), 0, 4e-3);
   unit_test_cmp('1d elec list equivalence (row/col)',unit_test_elec_pos(imdl1), unit_test_elec_pos(imdl2));
   unit_test_cmp('1d vs. 2d elec list equivalence',unit_test_elec_pos(imdl1), unit_test_elec_pos(imdl3));
   unit_test_cmp('2D PEM w/o node mashing, no vertical relief',[1:6]'*[1 0], unit_test_elec_pos(imdl3Hnm2d), eps);
   unit_test_cmp('2D PEM *is* PEM',length(imdl3Hnm2d.fwd_model.electrode(1).nodes),1)
   unit_test_cmp('3D PEM w/o node mashing, no vertical relief',[1:6]'*[1 0 0], unit_test_elec_pos(imdl3Hnm3d), eps);
   unit_test_cmp('3D PEM *is* PEM',length(imdl3Hnm3d.fwd_model.electrode(1).nodes),1)
   unit_test_cmp('1d vs. 2d + y=2 elec list equivalence',unit_test_elec_pos(imdl1), unit_test_elec_pos(imdl4)-([1:6]*0+2)'*[0 1]);
   unit_test_cmp('1d vs. 2d + y=2 - 135 deg elec eq', ...
                 unit_test_elec_pos(imdl1), ...
                 unit_test_elec_pos(imdl5, R(135), -X), eps*10);
   unit_test_cmp('2d with vertical geometry (mash nodes) CEM', ...
                 elec_pos_2d, ...
                 unit_test_elec_pos(imdl2dc), 0.01);
   unit_test_cmp('3d with vertical geometry (mash nodes) CEM', ...
                 elec_pos_3d, ...
                 unit_test_elec_pos(imdl3dc), 0.001);
   unit_test_cmp('2d with vertical geometry (mash nodes) PEM', ...
                 elec_pos_2d, ...
                 unit_test_elec_pos(imdl2dp), 10*eps);
   unit_test_cmp('3d with vertical geometry (mash nodes) PEM', ...
                 elec_pos_3d, ...
                 unit_test_elec_pos(imdl3dp), 10*eps);
   unit_test_cmp('2d with vertical geometry (w/o mash nodes) PEM', ...
                 elec_pos_2d, ...
                 unit_test_elec_pos(imdl2df), -10*eps);
   unit_test_cmp('3d with vertical geometry (w/o mash nodes) PEM', ...
                 elec_pos_3d, ...
                 unit_test_elec_pos(imdl3df), -10*eps);

   unit_test_cmp('h32a dual 3d/2d', ...
                 elec_pos_3d, ...
                 unit_test_elec_pos(imdlh32), 0.001);
   unit_test_cmp('H32a dual 3d/2d', ...
                 elec_pos_3d, ...
                 unit_test_elec_pos(imdlH32), 10*eps);
   unit_test_cmp('h22a dual 2d/2d', ...
                 elec_pos_2d, ...
                 unit_test_elec_pos(imdlh22), 0.002);
   unit_test_cmp('H22a dual 2d/2d', ...
                 elec_pos_2d, ...
                 unit_test_elec_pos(imdlH22), 10*eps);

clf; subplot(221); show_fem(imdl1.fwd_model); title('models match? A');
     subplot(222); show_fem(imdl5.fwd_model); title('models match? C');
     subplot(223); show_fem(imdl2dc.fwd_model); title('2d deformations');
     subplot(224); show_fem(imdl3dc.fwd_model); title('3d deformations'); view([0 -1 0.01]);

if 1
   clf; subplot(221); show_fem(imdlh22.fwd_model);
        subplot(222); show_fem(imdlh32.fwd_model);
        subplot(223); show_fem(imdlH22.fwd_model);
        subplot(224); show_fem(imdlH32.fwd_model);
end

if 0
   clf; subplot(221); show_fem(imdlh22.rec_model);
        subplot(222); show_fem(imdlh32.rec_model);
        subplot(223); show_fem(imdlH22.rec_model);
        subplot(224); show_fem(imdlH32.rec_model);
end

   imdl = mk_geophysics_model('h2p5a', ne, {'extra_ng_code', 'solid tt = orthobrick(-1,-1,-1;-0.5,0,-0.5);\ntlo tt;\n'});
   clf; show_fem(imdl.fwd_model);
   % Note: this just has to be swallowed safely...

function xyz = unit_test_elec_pos(imdl, R, X)
   if nargin < 2; R = 1; end
   if nargin < 3; X = 0; end
   fmdl = imdl.fwd_model;
   xyz = zeros(length(fmdl.electrode),size(fmdl.nodes,2))*NaN;
   for i = 1:length(fmdl.electrode)
      nn = fmdl.nodes(fmdl.electrode(i).nodes,:);
      xyz(i,:) = (max(nn,[],1) + min(nn,[],1))/2;
      xyz(i,:) = xyz(i,:)*R + X;
   end
