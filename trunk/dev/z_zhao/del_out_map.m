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
IN=bwperim(IN,8); %contour
xy = [X(IN)'; Y(IN)'];
out_ind=[];
max_x=max(xy(1,:));
max_y=max(xy(2,:));
min_x=min(xy(1,:));
min_y=min(xy(2,:));
% max_line_square=0;
% for i=1:length(xy(1,:))
%     temp=sum((repmat(xy(1:2,i),1,size(xy(:,:),2))-xy(:,:)).^2);
%     if any(temp>max_line_square)
%         max_line_square=max(temp);
%     end
% end

for i=1:length(xyzr)
% if round(xyzr(1,i)*10)==-7 && round(xyzr(2,i)*10)==4
%     M=M;
% end
    if xyzr(1,i)>max_x || xyzr(1,i)<min_x || xyzr(2,i)>max_y || xyzr(2,i)<min_y
        out_ind=[out_ind i];
        continue;
    end
    distance_to_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,:),2))-xy(:,:)).^2);
    on_boundary = distance_to_boundary <= xyzr(4,i)^2;
    if any(on_boundary)  
        out_ind=[out_ind i];
        if xyzr(1,i)<0 && xyzr(2,i)<0 % third
            ind=find(xyzr(1,:)<xyzr(1,i)&xyzr(2,:)<xyzr(2,i));
            out_ind=[out_ind ind];
        elseif xyzr(1,i)>0 && xyzr(2,i)<0 % fourth
            ind=find(xyzr(1,:)>xyzr(1,i)&xyzr(2,:)<xyzr(2,i));
            out_ind=[out_ind ind];
        elseif xyzr(1,i)>0 && xyzr(2,i)>0 % first
            ind=find(xyzr(1,:)>xyzr(1,i)&xyzr(2,:)>xyzr(2,i));
            out_ind=[out_ind ind];
        else  % second, xyzr(1,i)<0 && xyzr(2,i)>0
            ind=find(xyzr(1,:)<xyzr(1,i)&xyzr(2,:)>xyzr(2,i));
            out_ind=[out_ind ind];
        end
    end

%     out_boundary = distance_to_boundary >= max_line_square;
%     if any(out_boundary)  
%         out_ind=[out_ind i];
%     end
end
xyzr(:,out_ind)=[];
