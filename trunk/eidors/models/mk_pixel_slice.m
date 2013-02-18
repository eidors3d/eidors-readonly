function [imdl fmdl] = mk_pixel_slice(imdl,level,opt)
%MK_PIXEL_SLICE create a pixel model to reconstruct on
% OUT = MK_PIXEL_SLICE(MDL, LEVEL, OPT) creates a slice of pixels as a
% model to reconstruct on. 
%
% Inputs:
%  MDL   = an EIDORS forward or inverse model structure
%  LEVEL = either a vector of x-,y-, and z-intercepts of the cut plane or
%          a single value interpreted as the height of a single cut in the
%          z plane (by default, a horizontal slice at the average electrode
%          height will be created)
%  OPT   = an option structure with the following fields and defaults:
%     opt.imgsz = [32 32];    % dimensions of the pixel grid
%     opt.square_pixels = 0;  % adjust imgsz to get square pixels
%     opt.do_coarse2fine = 1; % calcuate c2f on the forward model
%     opt.z_depth = inf;      % z_depth to use with mk_coarse_fine_mapping
%
% Output depends on the type of model suplied. If MDL is a fwd_model
% structure then OUT is a rec_model. If MDL is an inv_model, then OUT is a
% modified version of it, with the pixel slice in inv_model.rec_model and
% updated inv_model.fwd_model.coarse2fine
% 
% [OUT FMDL] = MK_PIXEL_SLICE(MDL, ...) also returns the forward model
% structure with the coarse2fine field.
%
% See also MK_COARSE_FINE_MAPPING, MK_GRID_MODEL

% (C) 2013 Bartlomiej Grychtol. License: GPL version 2 or 3
% $Id$

if isstr(imdl) && strcmp(imdl,'UNIT_TEST'),do_unit_test;return;end;

switch(imdl.type)
    case 'inv_model'
        fmdl = imdl.fwd_model;
    case 'fwd_model'
        fmdl = imdl;
    otherwise
        error('An EIDORS inverse or forward model struct required');
end

if nargin < 2, opt = struct; end
if nargin > 1
   opt.level = level;
end
opt = parse_opts(fmdl,opt);

tmp = fmdl;
[NODES R T]=level_model(fmdl,opt.level);
tmp.nodes = NODES';
slc = mdl_slice_mesher(tmp,[inf inf 0]);
slc = slc.fwd_model;
mingrid = min(slc.nodes);
maxgrid = max(slc.nodes);
bnd = find_boundary(slc);
contour_boundary = order_loop(slc.nodes(unique(bnd),:));

if opt.square_pixels ==1
    mdl_sz = maxgrid - mingrid;
    mdl_AR = mdl_sz(1)/mdl_sz(2);
    img_AR = opt.imgsz(1)/opt.imgsz(2);
    if mdl_AR < img_AR
        delta = (mdl_sz(2) * img_AR - mdl_sz(1)) /2;
        mingrid(1) = mingrid(1) - delta;
        maxgrid(1) = maxgrid(1) + delta;
    elseif mdl_AR > img_AR
        delta = (mdl_sz(1)/img_AR - mdl_sz(2)) / 2;
        mingrid(2) = mingrid(2) - delta;
        maxgrid(2) = maxgrid(2) + delta;
    end
end

xgrid = linspace(mingrid(1),maxgrid(1),opt.imgsz(1)+1);
ygrid = linspace(mingrid(2),maxgrid(2),opt.imgsz(2)+1);
rmdl = mk_grid_model([],xgrid,ygrid);
x_pts = xgrid(1:end-1) + 0.5*diff(xgrid);
y_pts = ygrid(1:end-1) + 0.5*diff(ygrid);
% y_pts = fliplr(y_pts); %medical

% NOTE: This controls the image resolution. If you want higher res, you
% need to either specify it in opt.imgsz or manually overwrite (or remove)
% the imdl.rec_model.mdl_slice_mapper.
rmdl.mdl_slice_mapper.x_pts = x_pts;
rmdl.mdl_slice_mapper.y_pts = y_pts;
rmdl.mdl_slice_mapper.level = opt.level;
rmdl.mdl_slice_mapper.model_2d = 1;
x_avg = conv2(xgrid, [1,1]/2,'valid');
y_avg = conv2(ygrid, [1,1]/2,'valid');
[x,y] = ndgrid( x_avg, y_avg);
inside = inpolygon(x(:),y(:),contour_boundary(:,1),contour_boundary(:,2) );
ff = find(~inside);
rmdl.elems([2*ff, 2*ff-1],:)= [];
rmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
rmdl.coarse2fine(:,ff)= [];
% rmdl.boundary = find_boundary(rmdl);
% show individual elements (more like how the 2d grid models display)
rmdl.boundary = rmdl.elems;

if opt.do_coarse2fine
    % to calculate c2f, models must be aligned
    tmp = rmdl;
    tmp.mk_coarse_fine_mapping.f2c_offset  = T;
    tmp.mk_coarse_fine_mapping.f2c_project = R;
    tmp.mk_coarse_fine_mapping.z_depth     = opt.z_depth;
    fmdl.coarse2fine = mk_coarse_fine_mapping(fmdl,tmp);
end

rmdl.nodes(:,3) = 0;
rmdl.nodes =  (R\rmdl.nodes' + T'*ones(1,length(rmdl.nodes)))';
slc.nodes =  (R\slc.nodes' + T'*ones(1,length(slc.nodes)))';
if isfield(slc, 'electrode')
    for i = flipud(1:numel(slc.electrode))
        x_elec = slc.nodes( [slc.electrode(i).nodes], 1);
        y_elec = slc.nodes( [slc.electrode(i).nodes], 2);
        z_elec = slc.nodes( [slc.electrode(i).nodes], 3);
        elec(i).z_contact = slc.electrode(i).z_contact;
        elec(i).pos       = [x_elec, y_elec, z_elec];
    end
    rmdl.electrode = elec;
end

rmdl.show_slices.levels = opt.level;

switch imdl.type
   case 'inv_model'
      imdl.rec_model = rmdl;
      imdl.fwd_model = fmdl;
   case 'fwd_model'
      imdl = rmdl;
end
      
 
 function opt = parse_opts(fmdl, opt)
    if ~isfield(opt, 'imgsz'),     
        opt.imgsz = [32 32]; 
    end
    if ~isfield(opt, 'square_pixels')
        opt.square_pixels = 0;
    end
    if ~isfield(opt, 'level')
        opt.level = get_elec_level(fmdl);
    else
        if numel(opt.level) ==1
            opt.level = [inf inf opt.level];
        end
    end
    if ~isfield(opt, 'do_coarse2fine')
        opt.do_coarse2fine = 1;
    end
    if ~isfield(opt, 'z_depth')
        opt.z_depth = inf;
    end
    
 function [NODE R T] = level_model( fwd_model, level )

   vtx= fwd_model.nodes;
   [nn, dims] = size(vtx);
   if dims ==2 % 2D case
       NODE= vtx';
       return;
   end

   % Infinities tend to cause issues -> replace with realmax
   % Don't need to worry about the sign of the inf
   level( isinf(level) | isnan(level) ) = realmax;
   level( level==0 ) =     1e-10; %eps;

   % Step 1: Choose a centre point in the plane
   %  Weight the point by it's inv axis coords
   invlev= 1./level;
   ctr= invlev / sum( invlev.^2 );

   % Step 2: Choose basis vectors in the plane
   %  First is the axis furthest from ctr
   [jnk, s_ax]= sort( - abs(level - ctr) );
   v1= [0,0,0]; v1(s_ax(1))= level(s_ax(1));
   v1= v1 - ctr;
   v1= v1 / norm(v1);

   % Step 3: Get off-plane vector, by cross product
   v2= [0,0,0]; v2(s_ax(2))= level(s_ax(2));
   v2= v2 - ctr;
   v2= v2 / norm(v2);
   v3= cross(v1,v2);

   % Step 4: Get orthonormal basis. Replace v2
   v2= cross(v1,v3);

   % Step 5: Get bases to point in 'positive directions'
   v1= v1 * (1-2*(sum(v1)<0));
   v2= v2 * (1-2*(sum(v2)<0));
   v3= v3 * (1-2*(sum(v3)<0));
   
   R = [v1;v2;v3];
   T = ctr;

   NODE= R * (vtx' - T'*ones(1,nn) );

function elec_lev = get_elec_level(fmdl)
    z_elec= fmdl.nodes( [fmdl.electrode(:).nodes], 3);
    min_e = min(z_elec); max_e = max(z_elec);
    elec_lev = [inf,inf,mean([min_e,max_e])];
    
function do_unit_test
    imdl = mk_common_model('n3r2',[16,2]);
    opt.square_pixels = 1;
    opt.imgsz = [16 16];
    mdl = mk_pixel_slice(imdl.fwd_model,[inf 2 2.5], opt);
    img = mk_image(mdl,1);
    
    subplot(221)
    show_fem(imdl.fwd_model);
    view([-50 10])

    subplot(222)
    show_fem(img);
    zlim([0 3]);
    ylim([-1 1])
    xlim([-1 1]);
    view([-50 10])
    
    subplot(223)
    show_slices(img);
    
    subplot(224)
    imdl = mk_pixel_slice(imdl);
    img = mk_image(imdl.rec_model,1);
    show_fem(img);
    zlim([0 3]);
    ylim([-1 1])
    xlim([-1 1]);
    view([-50 10])
    