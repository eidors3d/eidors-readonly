function PSF = new_desired_soln(xyz,radius, opt)


if ischar(xyz) && strcmp(xyz,'UNIT_TEST'), do_unit_test; return, end

[xyzr, radius, opt] = parse_opt(xyz, radius, opt);

mdl = opt.rec_model;
show_fem(mdl)

copt.cache_obj = {opt};
copt.fstr = 'new_desired_soln';
PSF = eidors_cache(@do_new_desired_soln,{xyzr, radius, opt},copt);


end


function PSF = do_new_desired_soln(xyz, radius, opt)
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
    th = log(1e4)/opt.steepness;
    for i=1:size(xyz,2);
%         if i == 101, error('jnk'); end
        progress_msg(i,num_it);
        farel = far_elems(Xnodes,Ynodes,Znodes, xyz(:,i), radius(i), th);
        el_idx = find(~farel);
        X = []; Y = []; Z = [];
        idx = []; factor = [];
        for e = el_idx'
            [x,y,z] = interp_elem(mdl,e,radius(i), opt);
            X = [X ; x(:)]; Y = [Y ; y(:)]; Z = [Z ; z(:)];
            n = numel(x);
            idx = [idx; e*ones(n,1)];
            factor = [factor; ones(n,1)/n];
        end
        D = sqrt(sum(bsxfun(@minus,[X Y Z],xyz(:,i)').^2, 2));
        x = D - radius(i);
        tmp = 1 ./ (1 + exp( opt.steepness * x));
        PSF(:,i) = full(sparse(idx,1,tmp(:) .* factor,size(mdl.vox,1),1));
    end
    progress_msg(Inf);
end

function [X, Y, Z] = voxnodes(mdl)
    x = mdl.nodes(:,1); X = x(mdl.vox);
    y = mdl.nodes(:,2); Y = y(mdl.vox);
    try
        z = mdl.nodes(:,3); Z = z(mdl.vox);
    catch
        Z = [];
    end
        
    
end


function [x,y,z,n] = interp_elem(mdl,e,radius, opt)
    maxsep = radius/5;

    minnode = min(mdl.nodes(mdl.vox(e,:),:));
    maxnode = max(mdl.nodes(mdl.vox(e,:),:));
    vec = {};
    for d = 1:opt.n_dim
        sep = maxnode(d) - minnode(d);
        N = max(3, ceil(sep/maxsep)+1);
        v = linspace(minnode(d),maxnode(d),N);
        v(end) = [];
        v = v + .5*(maxnode(d)-minnode(d))/(N-1);
        vec{d} = v;
    end
    switch opt.n_dim
        case 3
            [x, y, z] = ndgrid(vec{:});
        case 2
            [x, y] = ndgrid(vec{:});
            z = [];
    end
%     clf
%     show_fem(mdl);
%     hold on
%     plot3(x(:),y(:),z(:),'.');
%     hold off
%     pause
end

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


function [xyz, radius, opt] = parse_opt(xyz, radius, opt)

    scale_radius = false;
    if isempty(radius)
        radius = xyz(end,:);
        scale_radius = true;
        xyz(end,:) = [];
    end
    mdl = opt.rec_model; % must exist
    opt.n_dim = size(mdl.nodes,2);
    xyz = xyz(1:opt.n_dim, :); % ingore z if model is 2D
    

    opt.meshsz = [];
    try
        for i = 1:3
            opt.meshsz = [opt.meshsz min(mdl.nodes(:,i)) max(mdl.nodes(:,i))];
        end
    end
    
    opt.n_dim = length(opt.meshsz)/2;
    opt.meshsz = reshape(opt.meshsz,2,[])';
%     if size(opt.meshsz, 1)==2
%         opt.meshsz(3,:) = 0;
%     end
    opt.minpt = opt.meshsz(:,1);
    opt.maxpt = opt.meshsz(:,2);
    opt.range = opt.maxpt - opt.minpt;
    opt.maxrange = max(opt.range(1:2));

    %rescale points to between -1 and 1 in the x-y plane
    xyz = 2*bsxfun(@minus, xyz, mean(opt.meshsz,2))/opt.maxrange;
    mdl.nodes = 2*bsxfun(@minus, mdl.nodes, mean(opt.meshsz(1:size(mdl.nodes,2),:),2)')/opt.maxrange;
    opt.rec_model = mdl;
    if scale_radius
        radius = radius / opt.maxrange;
    end
    
    if ~isfield(opt, 'steepness')
        opt.steepness = 10/mean(radius);
    end
    
    if numel(radius) == 1
        radius = ones(1,size(xyz,2)) * radius;
    end
end

function do_unit_test
    v = linspace(-1,1,32);
    mdl= mk_grid_model([],v,v,[0 .7 1:.2:2 2.3 3]);
    opt.rec_model = mdl;
    xyzr = zeros(5,4);
    xyzr(:,4) = .25;
    xyzr(:,3) = .5:.5:2.5;
    sol = new_desired_soln(xyzr',[],opt);
    img = mk_image(mdl,0);
    for i = 1:5
        subplot(2,3,i)
        img.elem_data= sol(:,i);
        show_3d_slices(img,xyzr(i,3),xyzr(i,2),xyzr(i,1));
    end
end