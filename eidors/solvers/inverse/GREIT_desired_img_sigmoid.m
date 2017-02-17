function PSF = GREIT_desired_img_sigmoid(xyz,radius, opt)
%GREIT_DESIRED_IMG_SIGMOID sigmoid-decay desired image function for GREIT 
% PSF= GREIT_desired_img_sigmoid(xyc, radius, opt)
%   xyz     - array of centers of desired target images
%               [2xN] - xy only
%               [3xN] - xyz or, if radius = [], xyr
%               [4xN] - xyzr (radius is ignored)
%   radius  - the radius of the target on the desired image as a fraction
%             of the model radius (half the larger dimension in xy)
%   opt     - a struct with these fields:
%      .rec_model   a 2D or 3D model as generated by MK_GRID_MODEL that may
%                   have had some pixels/voxels removed (a rec_model that 
%                   is not a complete rectangle/parallelepiped). 
%      .steepness   [optional] a factor controling the amount of blur, see
%                   below for details. Lower steepness gives more blur, but
%                   if the value is too low, image may not reach the value 
%                   of 1 at the center of the target. 
%                   May be specified as a scalar, a [1xN] vector, or a
%                   function handle with the signature:
%                        func(pts)
%                   where pts is either [2xN] or [3xN], dependig on the xyz
%                   function input, e.g. @(xyz) 50./(xyz(3,:))
%                   Default: 10./radius
%      .desired_img_radius
%                   [optional] Overwrites the radius input. May be 
%                   specified as a scalar, a [1xN] vector or a function 
%                   handle with the signature:
%                        func(pts)
%                   where pts is either [2xN] or [3xN], dependig on the xyz
%                   function input, e.g. @(xyz) abs(xyz(3,:))/5
%      .threshold   [optional] voxels where the function is smaller than
%                   threshold will be set 0. Values bigger than 1-threshold
%                   will be 1. The smaller the threshold the more
%                   computationally expensive is the evaluation. 
%                   Default: 1e-4;
%
% The desired images approximate in each pixel the area integral of:
%       f(r) = 1 / (1 + exp(s*(|r-r0| - radius)))
% where
%       r   - position vector in 2D/3D space
%       s   - opt.steepness
%       r0  - target center
% For |r-r0| = radius, f(r) = 0.5.
% 
% As of 2015-03-29, this is the default desired image function used by
% MK_GREIT_MODEL.
%
% See also: CALC_GREIT_RM, MK_GREIT_MODEL, MK_PIXEL_SLICE

% (C) 2015 Bartlomiej Grychtol - all rights reserved Swisstom AG.
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<

if ischar(xyz) && strcmp(xyz,'UNIT_TEST'), do_unit_test; return, end

[xyzr, radius, opt] = parse_opt(xyz, radius, opt);

copt.cache_obj = {xyzr, radius, opt.rec_model.nodes, opt.rec_model.elems, opt.steepness};
copt.fstr = 'GREIT_desired_img_sigmoid';
PSF = eidors_cache(@desired_soln,{xyzr, radius, opt},copt);


end

%-------------------------------------------------------------------------%
% The main function
function PSF = desired_soln(xyz, radius, opt)
   num_it = size(xyz,2);
    progress_msg('Composing desired images:',0,num_it);
    mdl = opt.rec_model;
    switch opt.n_dim
        case 3
            mdl.vox = [mdl.elems(1:6:end,:) mdl.elems(6:6:end,:)];
        case 2
            mdl.vox = [mdl.elems(1:2:end,:) mdl.elems(2:2:end,:)];
    end
    [Xnodes,Ynodes,Znodes] = voxnodes(mdl);
    
    min_vox_edge = min( [min(diff(unique(mdl.nodes(:,1)))), ...
                        min(diff(unique(mdl.nodes(:,2)))), ...
                        min(diff(unique(mdl.nodes(:,3))))] );

    
    warned = false;
    interp_elem_new('reset',size(mdl.vox,1));
    
    % estimate the amount of memory needed to store PSF
    n_el = size(mdl.vox,1);
    n_pt = size(xyz,2);
    n_ep = n_el*n_pt;
    ll = eidors_msg('log_level',0);
    vox_box = sum(get_elem_volume(mdl));
    eidors_msg('log_level',ll);
    tgt_box = (2*mean(radius+opt.threshold./opt.steepness))^3;
    n_nz = ceil(n_ep*tgt_box/vox_box);
    SPARSE_STEP = ceil(n_nz/5);
    PSF = spalloc(n_el,n_pt,n_nz);
    
    for i=1:size(xyz,2)
        if nnz(PSF) > n_nz
           n_nz = nnz(PSF) + SPARSE_STEP;
           [r,c,v] = find(PSF);
           PSF = sparse(r,c,v,n_el,n_pt,n_nz); 
        end
        th = opt.threshold/opt.steepness(i);
        progress_msg(i,num_it);
        farel = far_elems(Xnodes,Ynodes,Znodes, xyz(:,i), radius(i), th);
        % also check for elements so close to the center they are approx 1
        close_el = false(size(farel));
        if radius(i) - th > min_vox_edge
           close_el = close_elems(Xnodes,Ynodes,Znodes, xyz(:,i), radius(i), th);
        end
        el_idx = find(~farel & ~close_el);
        % start with an initial size
        STEP = 1000;
        idx = zeros(STEP,1);
        factor = zeros(STEP,1);
        X = zeros(STEP,1);
        Y = zeros(STEP,1);
        Z = zeros(STEP,1);
        L = STEP;
        N = 1;

        for e = el_idx'
            [x,y,z] = interp_elem_new(mdl,e,radius(i), opt);
            n = numel(x);
            N1 = N; 
            N = N+n;
            if N > L % expand arrays as needed
                fill = zeros(STEP,1);
                idx = [idx; fill];
                factor = [factor; fill];
                X = [X; fill]; Y = [Y; fill]; Z = [Z; fill];
                L = L + STEP;
            end
            v = N1:N-1;
            idx(v) = e;
            factor(v) = 1/n;
            X(v) = x; Y(v) = y; Z(v) = z;
        end
        % remove spare length
        v = N:L;
        idx(v) = [];
        factor(v) = [];
        X(v) = []; Y(v) = []; Z(v) = [];
        
        if isempty(X)
            if ~warned
               warning('EIDORS:OutsidePoint',...
                  'Desired image generation failed for point %d (and maybe others)',i);
               warned = true;
            end
            % PSF(:,i) = 0; % not needed
            continue;
        end
        V = [X Y Z];
        D = sqrt(sum(bsxfun(@minus,V(:,1:opt.n_dim),xyz(1:opt.n_dim,i)').^2, 2));
        x = D - radius(i);
        tmp = 1 ./ (1 + exp( opt.steepness(i) * x));
        PSF(:,i) = sparse(idx,1,tmp(:) .* factor,size(mdl.vox,1),1);
        PSF(close_el,i) = 1;
    end
    progress_msg(Inf);
end

%-------------------------------------------------------------------------%
% Generate internal points in elements
function [X, Y, Z] = voxnodes(mdl)
    x = mdl.nodes(:,1); X = x(mdl.vox);
    y = mdl.nodes(:,2); Y = y(mdl.vox);
    try
        z = mdl.nodes(:,3); Z = z(mdl.vox);
    catch
        Z = [];
    end  
end

function [x,y,z] = interp_elem_new(mdl,e,radius,opt)
    persistent N_entries X  Y Z MAP N minnode maxnode done_elems sep
    if ischar(mdl) && strcmp(mdl,'reset')
        N_entries = 0; X = []; Y = []; Z = []; MAP = [];
        N       = ones(1,3);
        minnode = zeros(e,3);
        maxnode = zeros(e,3);
        sep     = zeros(e,3);
        done_elems = false(e,1);
        return;
    end
    maxsep = radius/5;
    
    if ~done_elems(e)
       minnode(e,:) = min(mdl.nodes(mdl.vox(e,:),:));
       maxnode(e,:) = max(mdl.nodes(mdl.vox(e,:),:));
       sep(e,:) = maxnode(e,:) - minnode(e,:);
       done_elems(e) = true;
    end
    
    N(1:opt.n_dim) = max(3, ceil(sep(e,:)/maxsep)+1);
    
    try
        entry = MAP(N(1),N(2),N(3));
        x = X{entry};
        y = Y{entry};
        z = Z{entry};
    catch
        entry = N_entries+1;
        
        vx = linvec(0,1,N(1));
        vx = vx + .5*(vx(2) -vx(1));
        
        vy = linvec(0,1,N(2));
        vy = vy + .5*(vy(2) -vy(1));
        
        
        switch opt.n_dim
            case 2
                [x, y] = grid2d(vx,vy);
                z = zeros(size(x));
            case 3
                vz = linvec(0,1,N(3));
                vz = vz + .5*(vz(2) -vz(1));
                [x, y, z] = grid3d(vx,vy,vz);
        end
        X{entry} = x;
        Y{entry} = y;
        Z{entry} = z;
        N_entries = entry;
        MAP(N(1),N(2),N(3)) = entry;
    end
    x = x*sep(e,1) + minnode(e,1);
    y = y*sep(e,2) + minnode(e,2);
    if opt.n_dim == 3
       z = z*sep(e,3) + minnode(e,3);
    end
    
end

%-------------------------------------------------------------------------%
% Generate internal points in elements
function [x,y,z] = interp_elem(mdl,e,radius, opt)
    maxsep = radius/5;

    minnode = min(mdl.nodes(mdl.vox(e,:),:));
    maxnode = max(mdl.nodes(mdl.vox(e,:),:));
    
    sep = maxnode - minnode;
    N = max(3, ceil(sep/maxsep)+1);

    vx = linvec(minnode(1),maxnode(1),N(1));
    vx = vx + .5*(vx(2) -vx(1));
        
    vy = linvec(minnode(2),maxnode(2),N(2));
    vy = vy + .5*(vy(2) -vy(1));
    

    switch opt.n_dim
        case 2
            [x, y] = grid2d(vx,vy);
            z = [];
        case 3
            vz = linvec(minnode(3),maxnode(3),N(3));
            vz = vz + .5*(vz(2) -vz(1));
            [x, y, z] = grid3d(vx,vy,vz);
    end
end


function out = linvec(v1, v2, N)
% equivalent to out= linspace(v1,v2,N); out(end) = []; minus error checking
out = v1 + (0:(N-2)).*(v2-v1)/(N-1);

end

function [x, y] = grid2d(vx, vy)
% ndgrid minus overhead, asumes both vx and vy are row vectors
    vx = vx.';
    x = repmat(vx,size(vy));
    y = repmat(vy,size(vx));
end

function [x, y, z] = grid3d(vx,vy,vz)
% ndgrid minus overhead
    sz = [numel(vx) numel(vy) numel(vz)];
    x = reshape(vx,[sz(1) 1 1]);
    x = repmat(x,[1 sz(2) sz(3)]);
    
    y = reshape(vy,[1 sz(2) 1]);
    y = repmat(y,[sz(1) 1 sz(3)]);
    
    z = reshape(vz,[1 1 sz(3)]);
    z = repmat(z,[sz(1) sz(2) 1]);
    

end


%-------------------------------------------------------------------------%
% Find elements where the function value is negligable
function farel = far_elems(Xnodes,Ynodes,Znodes,xyz,radius, th)
    farel = false(size(Xnodes,1),1);

    nodes_test = Xnodes < xyz(1) - radius - th;
    farel = farel | all(nodes_test,2);
    if all(farel), return, end;
    nodes_test = Xnodes > xyz(1) + radius + th;
    farel = farel | all(nodes_test,2);
    if all(farel), return, end;
    nodes_test = Ynodes < xyz(2) - radius - th;
    farel = farel | all(nodes_test,2);
    if all(farel), return, end;
    nodes_test = Ynodes > xyz(2) + radius + th;
    farel = farel | all(nodes_test,2);
    if all(farel), return, end;
    if ~isempty(Znodes)
        nodes_test  = Znodes > xyz(3) + radius + th;
        farel = farel | all(nodes_test,2);
        if all(farel), return, end;
        nodes_test  = Znodes < xyz(3) - radius - th;
        farel = farel | all(nodes_test,2);
        if all(farel), return, end;
    end
    idx = find(~farel);
end

%-------------------------------------------------------------------------%
% Find elements where the function value is essentially 1
function farel = close_elems(Xnodes,Ynodes,Znodes,xyz,radius, th)
  
   farel = true(size(Xnodes,1),1);

    nodes_test = Xnodes < xyz(1) + radius - th;
    farel = farel & all(nodes_test,2);
    if ~any(farel), return, end;
    nodes_test = Xnodes > xyz(1) - radius + th;
    farel = farel & all(nodes_test,2);
    if ~any(farel), return, end;
    nodes_test = Ynodes < xyz(2) + radius - th;
    farel = farel & all(nodes_test,2);
    if ~any(farel), return, end;
    nodes_test = Ynodes > xyz(2) - radius + th;
    farel = farel & all(nodes_test,2);
    if ~any(farel), return, end;
    if ~isempty(Znodes)
        nodes_test  = Znodes > xyz(3) - radius + th;
        farel = farel & all(nodes_test,2);
        if ~any(farel), return, end;
        nodes_test  = Znodes < xyz(3) + radius - th;
        farel = farel & all(nodes_test,2);
        if ~any(farel), return, end;
    end
    idx = find(farel);
end


%-------------------------------------------------------------------------%
% Parse options
function [xyz, radius, opt] = parse_opt(xyz, radius, opt)

    scale_radius = false;
    if isempty(radius)
        radius = xyz(end,:);
        scale_radius = true;
        xyz(end,:) = [];
    end
    
    if isfield(opt,'desired_img_radius')
       scale_radius = false;
       if isnumeric(opt.desired_img_radius)
          radius = opt.desired_img_radius;
       end
       if isa(opt.desired_img_radius, 'function_handle')
          radius = feval(opt.desired_img_radius,xyz);
       end
    end
        
 
    
    mdl = opt.rec_model; % must exist
    opt.n_dim =mdl_dim(mdl);
    xyz = xyz(1:opt.n_dim, :); % ingore z if model is 2D

    

    opt.meshsz = [];
    try
        for i = 1:3
            opt.meshsz = [opt.meshsz min(mdl.nodes(:,i)) max(mdl.nodes(:,i))];
        end
    end
    
    opt.n_dim = length(opt.meshsz)/2;
    opt.meshsz = reshape(opt.meshsz,2,[])';

    opt.minpt = opt.meshsz(:,1);
    opt.maxpt = opt.meshsz(:,2);
    opt.range = opt.maxpt - opt.minpt;
    opt.maxrange = max(opt.range(1:2));

    %rescale points to between -1 and 1 in the x-y plane
    xyz = 2*bsxfun(@minus, xyz, mean(opt.meshsz,2))/opt.maxrange;
    mdl.nodes = 2*bsxfun(@minus, mdl.nodes, mean(opt.meshsz(1:size(mdl.nodes,2),:),2)')/opt.maxrange;
    opt.rec_model = mdl;
    if scale_radius
        radius = 2*radius / opt.maxrange;
    end

    if ~isfield(opt, 'steepness')
        opt.steepness = 10./radius;
    elseif isa(opt.steepness,'function_handle')
        opt.steepness = feval(opt.steepness,xyz);
    end
    
    if numel(opt.steepness) == 1
       opt.steepness = ones(1,size(xyz,2)) * opt.steepness;
    end
    
    if numel(radius) == 1
        radius = ones(1,size(xyz,2)) * radius;
    end
    
    if opt.n_dim == 2
       xyz(3,:) = 0;
    end
   
    if ~isfield(opt, 'threshold')
       opt.threshold = 1e-4;
    end
    
    opt.threshold = log(1/opt.threshold - 1); 
    
end

%-------------------------------------------------------------------------%
% UNIT_TEST
function do_unit_test
    eidors_cache off GREIT_desired_img_sigmoid
    v = linspace(-1,1,32);
    mdl= mk_grid_model([],v,v,[0 .7 1:.2:2 2.3 3]);
    opt.rec_model = mdl;
    opt.steepness = @(xyz) 50./xyz(3,:);
    opt.desired_img_radius = @(xyz) xyz(3,:)/5;
%     xyzr = zeros(5,4);
    xyzr(:,3) = .5:.5:2.5;
    xyzr(:,4) = .25;
    sol = GREIT_desired_img_sigmoid(xyzr',[],opt);
    img = mk_image(mdl,0);
    for i = 1:5
        subplot(2,3,i)
        img.elem_data= sol(:,i);
        show_3d_slices(img,xyzr(i,3),xyzr(i,2),xyzr(i,1));
    end
    eidors_msg('UNIT_TEST: Showed %d images to verify',i,0);
    eidors_cache on GREIT_desired_img_sigmoid
end
