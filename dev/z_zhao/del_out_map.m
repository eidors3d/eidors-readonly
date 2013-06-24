function [xyzr]=del_out_map(imgs,N,xyzr)
bnd_nodes = unique(imgs.fwd_model.boundary);
min_bb = min(imgs.fwd_model.nodes(bnd_nodes,:));
max_bb = max(imgs.fwd_model.nodes(bnd_nodes,:));
xspace = linspace(min_bb(1),max_bb(1),N);
yspace = linspace(min_bb(2),max_bb(2),N);
[X Y] = meshgrid(xspace,yspace);
imgs.calc_colours.npoints = N;
M = calc_slices(imgs,1);
IN = M==1;
IN=bwperim(IN,8); %contur
xy = [X(IN)'; Y(IN)'];
out_ind=[];
max_x=max(xy(1,:));
max_y=max(xy(2,:));
min_x=min(xy(1,:));
min_y=min(xy(2,:));
for i=1:length(xyzr)

    if xyzr(1,i)>max_x || xyzr(1,i)<min_x || xyzr(2,i)>max_y || xyzr(2,i)<min_y
        out_ind=[out_ind i];
    end
    
    out_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,:),2))-xy(:,:)).^2) <= xyzr(4,i)^2;
    if any(out_boundary)  
        out_ind=[out_ind i];
    end
    
end
xyzr(:,out_ind)=[];
