function h = show_3d_slices(img, varargin);
% show_3d_slices(img, z_cuts, x_cuts, y_cuts, any_cuts)
% Show a 3d view of an object with many slices through it
%  z_cuts = planes in z to do a cut
%  x_cuts = planes in x to do a cut
%  y_cuts = planes in y to do a cut
%  any_cuts = any single-slice specification accepted by LEVEL_MODEL_SLICE
% Default show 2 z_cuts and 1 x and 1 y cut
%
% h = show_3d_slices(...) returns a cell array of handles to the individual
% slices { z_slices, x_slices, y_slices, any_slices }
%
% See also: level_model_slice, mdl_slice_mesher

% (C) 2007-2012 Andy Adler & Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

% BUGS: for 2D images, large electrodes are not correctly shown (AA, jun 2019)

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

ok= 0;
try; if strcmp(img.type, 'image');
   ok = 1;
end; end
if ~ok; 
   error('EIDORS: show_3d_slices requires image as first parameter');
end
if size(img.fwd_model.elems,2) == 3; %show_3d_slices requires volume model
   img = fake_3d_model( img );
end

%multi-parametrization
img = data_mapper(img);

qfi = warning('query','EIDORS:FirstImageOnly');
% Get data and to c2f mapping if required
elem_data = get_img_data(img);
% check size
if size(elem_data,2) > 1
   q = warning('query','backtrace');
   warning('backtrace','off');
   warning('EIDORS:FirstImageOnly','show_3d_slices only shows first image');
   warning('backtrace',q.state);
   warning('off','EIDORS:FirstImageOnly');
end

% need to make sure boundary is just the outside
img.fwd_model.boundary = find_boundary(img.fwd_model);


[jnk,ref_lev,max_scale] = scale_for_display( elem_data);
try 
    img.calc_colours.ref_level; 
catch
    img.calc_colours.ref_level = ref_lev;
end
try
    img.calc_colours.clim;
catch
    img.calc_colours.clim = max_scale;
end
try 
    np = img.calc_colours.npoints;
catch
    np = calc_colours('npoints');
end
%  show_fem(img.fwd_model);
[x_cuts, y_cuts, z_cuts, any_cuts] = get_cuts(img,varargin{:});

hz = []; hy = []; hx = []; ha = [];
for i= 1:length(z_cuts)
    hz(i) = show_slice(img,[inf inf z_cuts(i)]);
    hold on
end

for i= 1:length(x_cuts)
    hx(i) = show_slice(img,[x_cuts(i) inf inf]);
    hold on
end

for i= 1:length(y_cuts)
    hy(i) = show_slice(img,[inf y_cuts(i) inf]);
    hold on
end

for i= 1:size(any_cuts,1)
    ha(i) = show_slice(img,any_cuts(i,:));
    hold on
end
hold off

if nargout > 0
   h = {hz,hx,hy,ha};
end

warning(qfi.state,'EIDORS:FirstImageOnly');

function h = show_slice(img, level)
    slc = mdl_slice_mesher(img,level);
    try slc.calc_colours = img.calc_colours; end
    h = show_fem(slc);

% create an extruded 3D model of height
function img = fake_3d_model( img );
  fmdl = img.fwd_model;
  fnodes = fmdl.nodes;
  felems = fmdl.elems;
  maxd = max( max(fnodes) - min(fnodes) );
  upd = 0.1*maxd; % up and down from zero
  Nn = num_nodes(fmdl);
  Ne = num_elems(fmdl);
  fmdn = fmdl;
  fmdn.nodes =[[fnodes,-upd+0*fnodes(:,1)];
               [fnodes,+upd+0*fnodes(:,1)]];
  for i=1:num_elecs(fmdl)
     eni = fmdn.electrode(i).nodes;
     fmdn.electrode(i).nodes = [eni(:);eni(:)+Nn]';
     if isempty(eni) && isfield(fmdn.electrode(i),'faces')
        eni = unique(fmdn.electrode(i).faces);
        fmdn.electrode(i).nodes = [eni(:);eni(:)+Nn]';
     end
  end
  if isfield(fmdn.electrode,'faces'); % Horrible we need to code it like this
     fmdn.electrode = rmfield(fmdn.electrode,'faces');
  end
  % Doesn't need to be a conformal model, don't worry about adjacency
  fmdn.elems =[[felems(:,[1,2,3]), felems(:,[1    ])+Nn];
               [felems(:,[  2,3]), felems(:,[1,2  ])+Nn];
               [felems(:,[    3]), felems(:,[1,2,3])+Nn]];
  ed= img.elem_data;
  if isfield(fmdl,'coarse2fine');
     c2f = fmdl.coarse2fine;
     fmdn.coarse2fine = [c2f;c2f;c2f]; 
     img.fwd_model = fmdn;
  else
     img = mk_image(fmdn,[ed;ed;ed]);
  end

    
function [x_cuts, y_cuts, z_cuts, any_cuts] =  get_cuts(img, varargin)
   mdl_max= max(img.fwd_model.nodes);
   mdl_min= min(img.fwd_model.nodes);
   any_cuts = [];
   if nargin==1;
      % Default show 2 z_cuts and 1 x and 1 y cut
       x_cuts= linspace(mdl_min(1), mdl_max(1), 3); x_cuts([1,3])=[];
       y_cuts= linspace(mdl_min(2), mdl_max(2), 3); y_cuts([1,3])=[];
       z_cuts= linspace(mdl_min(3), mdl_max(3), 4); z_cuts([1,4])=[];
   elseif nargin==2;
       z_cuts= varargin{1};
       x_cuts= [];
       y_cuts= [];
   elseif nargin==3;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= [];
   elseif nargin==4;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= varargin{3};
   elseif nargin==5;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= varargin{3};
       any_cuts= varargin{4};
   else 
       error('too many inputs');
   end
   


function do_unit_test
   imdl = mk_common_model('n3r2',[16,2]);
   img = mk_image(imdl,1); vh= fwd_solve(img);
   load datacom.mat A B;
   img.elem_data(A) = 1.2;
   img.elem_data(B) = 0.8;
   vi = fwd_solve(img);
   imgr = inv_solve(imdl, vh, vi);
   imgr.calc_colours.npoints =64;
   show_3d_slices(img);
   calc_colours('transparency_thresh', 0.25);
   show_3d_slices(imgr,[1.5,2],[],[]);

   subplot(4,4,1);  show_3d_slices(imgr,[1,2],0,0.5);
   subplot(4,4,2);  show_3d_slices(img ,[1,2],0,0.5);
   cuts = [inf, -2.5, 1.5; inf, 10, 1.5];
   subplot(4,4,3);  show_3d_slices(img ,[],[],[],cuts );
   
   imgr.calc_colours.transparency_thresh = -1;
   img.calc_colours.transparency_thresh = -1;
   subplot(4,4,4);  show_3d_slices(imgr,[1,2],0.3,0.7);
   subplot(4,4,5);  
   show_fem(img.fwd_model);
   hold on
   show_3d_slices(img ,[1,2],0.3,0.7);
   axis off

   subplot(4,4,6);  
   img.elem_data = img.elem_data*[0,1];
   img.get_img_data.frame_select = 2;
   show_3d_slices(img ,[1],[],[]);
   
   view(10,18);


   subplot(4,4,7);
   imgs = img;        ed = img.elem_data(:,1);
   imgs.elem_data= [1;1.2;0.8]*[1];
   imgs.fwd_model.coarse2fine = [ed == img.elem_data(1),
                                 ed == img.elem_data(2),
                                 ed == img.elem_data(3)];
   show_3d_slices(img,[1,2],0,0.7);
  
   subplot(4,4,8);
   show_3d_slices(img,[1],[],[],[1,1,1]);


   subplot(4,4,9);  
   vopt.imgsz = [10,10,10];
   % Build a model that has a c2f. This one isn't very good,
   % because of transpose, but that's OK for the test
   [img.fwd_model,tmp] = mk_voxel_volume(imdl.fwd_model,vopt);
   img.fwd_model.coarse2fine = (tmp.coarse2fine * ...
                      img.fwd_model.coarse2fine')';
   show_3d_slices(img,[1,2],0,0.5);

   subplot(4,4,10);
   fmdl = getfield(mk_common_model('a2C2',4),'fwd_model');
   fmdl.mat_idx{1} = 1:4;
   fmdl = mat_idx_to_electrode(fmdl,{1});
   img= mk_image(fmdl,1);
   show_3d_slices(img,[0])
   
   subplot(4,4,11);
   img.fwd_model.electrode(1).faces=[28,29;29,30;30,31];
   show_3d_slices(img,[0])

   subplot(4,4,12)
% Without C2f
   imdl = mk_common_model('a2c0',8); img= mk_image(imdl);
   img.elem_data([10,12]) = 1.1;
   show_3d_slices(img,[0],[0.2],[0]);


   subplot(4,4,13)
% With C2f
   ed = img.elem_data;
   img.fwd_model.coarse2fine = [ed==1, ed> 1];
   img.elem_data= [1;1.1];
   show_3d_slices(img,[0],[0.2],[0]);

   subplot(4,4,14)
% With internal electrode
   extra={'ball','solid ball = sphere(0.5,0.5,1;0.1);'};
   fmdl= ng_mk_cyl_models(2,[4,1],[0.2],extra); 
   fmdl.mat_idx_to_electrode.nodes_electrode = true;
   fmdl= mat_idx_to_electrode(fmdl,{2});
   show_3d_slices(mk_image(fmdl,1),[1.0],[0.5],[0.5]);
    
   subplot(4,4,15)
% With internal electrode
   extra={'ball','solid ball = sphere(0.5,0.5,1;0.1);'};
   fmdl= ng_mk_cyl_models(2,[4,1],[0.2],extra); 
   fmdl.mat_idx_to_electrode.nodes_electrode = false;
   fmdl= mat_idx_to_electrode(fmdl,{2});
   show_3d_slices(mk_image(fmdl,1),[1.0],[0.5],[0.5]);
