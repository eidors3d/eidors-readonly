function imdl = mk_rec_model(imdl,imgsz)

fmdl = imdl.fwd_model;

mingrid = min(fmdl.nodes);
maxgrid = max(fmdl.nodes);
% if opt.square_pixels ==1
    mdl_sz = maxgrid - mingrid; 
    mdl_AR = mdl_sz(1)/mdl_sz(2);
    img_AR = imgsz(1)/imgsz(2);
    if mdl_AR < img_AR
        delta = (mdl_sz(2) * img_AR - mdl_sz(1)) /2;
        mingrid(1) = mingrid(1) - delta;
        maxgrid(1) = maxgrid(1) + delta;
    elseif mdl_AR > img_AR
        delta = (mdl_sz(1)/img_AR - mdl_sz(2)) / 2;
        mingrid(2) = mingrid(2) - delta;
        maxgrid(2) = maxgrid(2) + delta;
    end       
% end

xgrid = linspace(mingrid(1),maxgrid(1),imgsz(1)+1);
ygrid = linspace(mingrid(2),maxgrid(2),imgsz(2)+1);
rmdl = mk_grid_model([],xgrid,ygrid);
x_pts = xgrid(1:end-1) + 0.5*diff(xgrid);
y_pts = ygrid(1:end-1) + 0.5*diff(ygrid); 
y_pts = fliplr(y_pts); %medical
% NOTE: This controls the image resolution. If you want higher res, you
% need to either specify it in opt.imgsz or manually overwrite (or remove)
% the imdl.rec_model.mdl_slice_mapper.
rmdl.mdl_slice_mapper.x_pts = x_pts;
rmdl.mdl_slice_mapper.y_pts = y_pts;
rmdl.mdl_slice_mapper.level = [inf inf 0];
x_avg = conv2(xgrid, [1,1]/2,'valid');
y_avg = conv2(ygrid, [1,1]/2,'valid');
[x,y] = ndgrid( x_avg, y_avg);

    z_elec= fmdl.nodes( [fmdl.electrode(:).nodes], 3);
    min_e = min(z_elec); max_e = max(z_elec);
    elec_lev = [inf,inf,mean([min_e,max_e])];
    fmdl.mdl_slice_mapper = rmdl.mdl_slice_mapper;
    fmdl.mdl_slice_mapper.level = elec_lev;

    % calculate mat_idx for the rec_model
     img = mk_image(fmdl,1);
     for i = 2:length(fmdl.mat_idx);
         img.elem_data(fmdl.mat_idx{i}) = i;
     end
     
    slice = calc_slices(img,elec_lev);
    slice = flipud(slice);
    slice = slice';
    inside = ~isnan(slice);
    ff= find(inside==0);
    

    
 
%     mat = reshape([slice(:)'; slice(:)'],1,[]);
% some assumptions for the time being:
% - first mat_idx is the background
% - each mat_idx describes one continues region
    mat = ones(size(rmdl.elems,1)/2,1);
    for i = 2:length(unique(slice(~isnan(slice))))
        % only keep the biggest region
        BW = slice == i;
        CC = bwconncomp(BW,4);
        n_obj = CC.NumObjects;
        l = cellfun('length',CC.PixelIdxList);
        [jnk, p] = max(l);
        mat(CC.PixelIdxList{p}) = i;
    end
    mat = reshape([mat'; mat'],1,[]);
    mat([2*ff, 2*ff-1])= [];
    rmdl.elems([2*ff, 2*ff-1],:)= [];
    rmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
    rmdl.coarse2fine(:,ff)= [];
    for i = 1:max(mat)
        rmdl.mat_idx(i) = {find(mat==i)'};
    end
 nodes = rmdl.nodes;
 for i = flipud(1:numel(fmdl.electrode))
    x_elec = mean(fmdl.nodes( [fmdl.electrode(i).nodes], 1));
    y_elec = mean(fmdl.nodes( [fmdl.electrode(i).nodes], 2));
    dist = (nodes(:,1)-x_elec).^2 + (nodes(:,2)-y_elec).^2;
    [jnk,e_node]= min(dist);
    elec(i).z_contact = 0.001;
    elec(i).nodes     = e_node;
 end
 rmdl.electrode = elec;

 imdl.rec_model = rmdl;
 imdl.fwd_model.coarse2fine = mk_coarse_fine_mapping(fmdl,rmdl);
 
