
%% 1. create the body
xy = [0 0;  1 0; 1 1; 0 1];
fmdl = ng_mk_2d_model({xy, 0.02}, [4 -.5],[0.05 200]);
show_fem(fmdl)
SSMM = [1 2 1 2;];
fmdl.stimulation = stim_meas_list(SSMM,4);

%% 2. create the object
obj_radius = 0.05;
obj = ng_mk_cyl_models([0, obj_radius, 0.01],[],[]);
bnd_nodes = unique(find_boundary(obj));
obj_contour = order_loop(obj.nodes(bnd_nodes,:));

show_fem(obj)
hold on
plot(obj_contour(:,1),obj_contour(:,2),'b')
hold off
%% 3. simulate movement
N = 4;
sqrtN = ceil(sqrt(N));
t = linspace(0,pi, N+1); t(end) = []; t = t(:);
xy = 0.5 + .25 * [sin(t) cos(t)];

obj_mv = obj;

for i = 1:size(xy,1)
    tmp = mdl;
    obj_mv.nodes = obj.nodes + xy(i,:);
    
    % nodes in the area to be replaced
    idx = sum((mdl.nodes - xy(i,:)).^2,2) <= (2 * obj_radius).^2;
    rm_nodes = find(idx);
    
    % faces to be removed
    rm_elems = find(sum(idx(mdl.elems),2) == 3); 
    tmp.elems(rm_elems,:) = [];
    
    % the newly created boundary
    tmp = fix_model(tmp, struct('edges',true));
    bnd_edges = sum(idx(tmp.edges),2) == 2;
    bnd_nodes = unique(tmp.edges(bnd_edges,:));
    bnd_contour = order_loop(tmp.nodes(bnd_nodes,:));
    
    
    % filler model 
    obj_mv_contour = flipud(obj_contour + xy(i,:));
    fill_mdl = gmsh_mk_2d_model({bnd_contour, obj_mv_contour, 1});
    
    % we can, but don't have to, remove the unused nodes
    usednodes = unique(tmp.elems(:));
    nidx = zeros(num_nodes(tmp),1);
    nidx(usednodes) = 1:length(usednodes);
    tmp.nodes(nidx==0,:) = [];
    tmp.elems = nidx(tmp.elems);
    
    % merge meshes
    mrg_mdl = merge_meshes(tmp, obj_mv, fill_mdl);
    
    mrg_img = mk_image(mrg_mdl,0);
    for k = 2:3
        mrg_img.elem_data(mrg_mdl.mat_idx{k}) = k;
    end
    
    subplot(sqrtN,sqrtN,i)
    show_fem(mrg_img)
    
    show_fem(tmp)
    hold on
    hh = show_fem(obj_mv); set(hh, 'edgecolor','b');
    hh = show_fem(fill_mdl); set(hh, 'edgecolor','c');
    plot(bnd_contour(:,1),bnd_contour(:,2),'r')
    hold off
    
    
end

