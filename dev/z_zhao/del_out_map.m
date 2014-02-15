function [xyzr]=del_out_map(imgs,N,xyzr)
% del_out_map: delete the similated targets out of the contour
%   [xyzr]=del_out_map(imgs,N,xyzr)
%
% Output: 
%   xyzr   - similated targets after correction
%
% Parameters:
%   imgs,    - img model on which to do simulations, 
%   N,       - reconstructed image size in pixels
%   xyzr     - similated targets

% Zhanqi Zhao 17.Sept.2013


bnd_nodes = unique(imgs.fwd_model.boundary);
min_bb = min(imgs.fwd_model.nodes(bnd_nodes,:));
max_bb = max(imgs.fwd_model.nodes(bnd_nodes,:));


% xspace = linspace(-max_xy,max_xy,N);
% yspace=xspace;
xspace = linspace(min_bb(1),max_bb(1),N);
% yspace = linspace(min_bb(2),max_bb(2),N); %64-64 whole
% yspace=linspace(-max_bb(1),-min_bb(1),N);
% yspace=linspace(min_bb(1),max_bb(1),N);
max_xy=max(abs([max_bb(1:2),min_bb(1:2)]));
yspace=linspace(-max_xy,max_xy/(max_bb(2)/-min_bb(2)),N);
[X Y] = meshgrid(xspace,-yspace);
imgs.calc_colours.npoints = N;
M = calc_slices(imgs,1);                    % only part, have to find a way
IN = M==1;
IN=bwperim(IN,8); %contour
xy = [X(IN)'; Y(IN)'];
% xy=xyzr(1:2,:);
out_ind=[];
max_x=max(xy(1,:));
max_y=max(xy(2,:));
min_x=min(xy(1,:));
min_y=min(xy(2,:));


for i=1:length(xyzr)
% if round(xyzr(1,i)*10)==-7 && round(xyzr(2,i)*10)==4
%     M=M;
% end

    if xyzr(1,i)>max_x || xyzr(1,i)<min_x || xyzr(2,i)>max_y || xyzr(2,i)<min_y
        out_ind=[out_ind i];
        continue;
    end
    distance_to_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,:),2))-xy(:,:)).^2);
%     on_boundary = distance_to_boundary <= xyzr(4,i)^2;
    on_boundary = distance_to_boundary <= (1.5*xyzr(4,i))^2; %inrease r for 50%
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

end
xyzr(:,out_ind)=[];
