function [imdl, fmdl] = mk_voxel_volume(varargin)
%MK_VOXEL_VOLUME create a voxel model to reconstruct on
% OUT = MK_VOXEL_VOLUME(MDL)
% OUT = MK_VOXEL_VOLUME(MDL, OPT)
%
% Inputs:
%  MDL   = an EIDORS forward or inverse model structure
%  OPT   = an option structure with the following fields and defaults:
%     opt.imgsz = [32 32 4]; % X, Y and Z dimensions of the voxel grid
%     opt.xvec  = []          % Specific X cut-planes between voxels
%                             % A scalar means the number of planes
%                             % Takes precedence over other options
%     opt.yvec  = []          % Specific Y cut-planes between voxels
%                             % A scalar means the number of planes
%                             % Takes precedence over other options
%     opt.zvec  = []          % Specific Z cut-planes between voxels
%                             % A scalar means the number of planes
%                             % Takes precedence over other options
%     opt.square_pixels = 0;  % adjust imgsz to get square pixels (in XY)
%     opt.cube_voxels = 0;    % adjust imgsz to get cube voxels (in XYZ)  
% NOT USED:
%     opt.do_coarse2fine = 1; % calcuate c2f on the forward model
%     opt.z_depth = inf;      % z_depth to use with mk_coarse_fine_mapping
%
% Output depends on the type of model suplied. If MDL is a fwd_model
% structure then OUT is a rec_model. If MDL is an inv_model, then OUT is a
% modified version of it, with the pixel slice in inv_model.rec_model and
% updated inv_model.fwd_model.coarse2fine
% 
% [OUT FMDL] = MK_VOXEL_VOLUME(MDL, ...) also returns the forward model
% structure with the coarse2fine field.
%
% See also MK_PIXEL_SLICE, MK_COARSE_FINE_MAPPING, MK_GRID_MODEL

% (C) 2015 Bartlomiej Grychtol - all rights reserved by Swisstom AG
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<

imdl = varargin{1};
% if input is 'UNIT_TEST', run tests
if ischar(imdl) && strcmp(imdl,'UNIT_TEST') 
    do_unit_test; clear imdl 
    return; 
end

if nargin < 2
    opt = struct;
else 
    opt = varargin{2};
end

switch(imdl.type)
    case 'inv_model'
        fmdl = imdl.fwd_model;
    case 'fwd_model'
        fmdl = imdl;
    otherwise
        error('An EIDORS inverse or forward model struct required');
end

opt = parse_opts(fmdl,opt);

copt.fstr = 'mk_voxel_volume';
copt.cache_obj = get_cache_obj(fmdl, opt);
[rmdl, c2f] = eidors_cache(@do_voxel_volume,{fmdl, opt},copt);

fmdl.coarse2fine = c2f;

switch imdl.type
   case 'inv_model'
      imdl.rec_model = rmdl;
      imdl.fwd_model = fmdl;
   case 'fwd_model'
      imdl = rmdl;
end


%-------------------------------------------------------------------------%
% The main function
function [rmdl, c2f] = do_voxel_volume(fmdl,opt)
%     opt.xvec = opt.xvec(11:12);
%     opt.yvec = opt.yvec(7:8);
%     opt.zvec = opt.zvec(1:2);
    rmdl = mk_grid_model([],opt.xvec,opt.yvec,opt.zvec);
%     fmdl.elems = fmdl.elems(6589,:);
    [c2f, m]  = mk_grid_c2f(fmdl, rmdl);
    inside = any(c2f);

    c2f(:,~inside) = [];
    rm = ~logical(rmdl.coarse2fine*inside');
    rmdl.elems(rm,:) = [];
    rmdl.coarse2fine(rm,:) = [];
    rmdl.coarse2fine(:,~inside) = [];
    
    
    bnd_fcs = ones(1,nnz(inside))*m.vox2face(inside,:) == 1;
    rmdl.boundary = m.faces(bnd_fcs,:);
    rmdl.inside   = inside; % the inside array is useful in other functions
    x_pts = opt.xvec(1:end-1) + 0.5*diff(opt.xvec);
    y_pts = opt.yvec(1:end-1) + 0.5*diff(opt.yvec);

    rmdl.mdl_slice_mapper.x_pts = x_pts;
    rmdl.mdl_slice_mapper.y_pts = y_pts;
    
    
%-------------------------------------------------------------------------%
% Assemble a reference object for caching
function cache_obj = get_cache_obj(fmdl, opt)
tmp = struct;
flds = {'nodes','elems'};
for f = flds;
    tmp.(f{1}) = fmdl.(f{1});
end
cache_obj = {tmp, opt};
    

%-------------------------------------------------------------------------%
% Parse option struct
 function opt = parse_opts(fmdl, opt)
    if ~isfield(opt, 'imgsz'),     
        opt.imgsz = [32 32 4]; 
    end
    if ~isfield(opt, 'square_pixels')
        opt.square_pixels = 0;
    end
    if ~isfield(opt, 'cube_voxels')
        opt.cube_voxels = 0;
    end
    if ~isfield(opt, 'xvec')
        opt.xvec = [];
    end
    if ~isfield(opt, 'yvec')
        opt.yvec = [];
    end
    if ~isfield(opt, 'zvec')
        opt.zvec = [];
    end
    if isempty(opt.xvec) && isempty(opt.imgsz)
        error('EIDORS:WrongInput', 'opt.imgsz must not be empty if opt.xvec is empty or absent');
    end
    if isempty(opt.yvec) && numel(opt.imgsz) < 2
        error('EIDORS:WrongInput', 'opt.imgsz must have at least 2 elements if opt.yvec is empty or absent');
    end
    if isempty(opt.zvec) && numel(opt.imgsz) < 3
        error('EIDORS:WrongInput', 'opt.imgsz must have 3 elements if opt.zvec is empty or absent');
    end
    
    mingrid = min(fmdl.nodes);
    maxgrid = max(fmdl.nodes);
    mdl_sz = maxgrid - mingrid;
    
    allempty = isempty(opt.xvec) & isempty(opt.yvec) & isempty(opt.zvec);
    if opt.cube_voxels && ~allempty
        warning('EIDORS:IncompatibleOptions','Option cube_voxels is ignored when xvec, yvec or zvec is specifed');
    end
    if opt.cube_voxels && allempty
        side_sz = max(mdl_sz ./ opt.imgsz(1:3));
        n_steps = ceil(mdl_sz / side_sz);
        mdl_ctr = mingrid + mdl_sz/2;
        mingrid = mdl_ctr - n_steps/2 * side_sz;
        maxgrid = mdl_ctr + n_steps/2 * side_sz;
        opt.imgsz = n_steps;
 
    elseif opt.square_pixels
        if ~isempty(opt.xvec) || ~isempty(opt.yvec)
            warning('EIDORS:IncompatibleOptions','Option square_pixels is ignored when xvec or yvec is specifed');
        else
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
    end
    if isempty(opt.xvec)
        opt.xvec = linspace(mingrid(1), maxgrid(1), opt.imgsz(1)+1);
    end      
    if isempty(opt.yvec)
        opt.yvec = linspace(mingrid(2), maxgrid(2), opt.imgsz(2)+1);
    end
    if isempty(opt.zvec)
        opt.zvec = linspace(mingrid(3), maxgrid(3), opt.imgsz(3)+1);
    end
    
    if ~isfield(opt, 'do_coarse2fine')
        opt.do_coarse2fine = 1;
    end
    if ~isfield(opt, 'z_depth')
        opt.z_depth = inf;
    end

%-------------------------------------------------------------------------%
% Perfom unit tests
function do_unit_test

fmdl = mk_library_model('adult_male_16el');
% fmdl= ng_mk_cyl_models([2,2,.4],[16,1],[.1,0,.025]);
opt.square_pixels = 1;
[rmdl, fmdl] = mk_voxel_volume(fmdl, opt);


subplot(121)
rimg = mk_image(rmdl,0);
rimg.elem_data = zeros(size(rmdl.coarse2fine,2),1);
idx = round(rand(5,1)* length(rimg.elem_data));
rimg.elem_data(idx) = 1;
show_fem(rimg);

subplot(122)
img = mk_image(fmdl,0);
img.elem_data = fmdl.coarse2fine * rimg.elem_data;
img.calc_colours.ref_level = 0;
img.calc_colours.transparency_thresh = 1e-2;

show_fem(img);





