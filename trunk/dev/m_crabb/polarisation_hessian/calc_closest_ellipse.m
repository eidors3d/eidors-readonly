function [ fmdl ] = calc_closest_ellipse( fmdl )
%CALC_CLOSEST_ELLIPSE 
%   Heuristic closest ellipse, based on longest axis of triangle elements

%
lens1 = sqrt(sum((fmdl.nodes(fmdl.elems(:,2),:) - fmdl.nodes(fmdl.elems(:,1),:)).^2,2));
lens2 = sqrt(sum((fmdl.nodes(fmdl.elems(:,3),:) - fmdl.nodes(fmdl.elems(:,1),:)).^2,2));
lens3 = sqrt(sum((fmdl.nodes(fmdl.elems(:,3),:) - fmdl.nodes(fmdl.elems(:,2),:)).^2,2));

%
len = [lens1, lens2, lens3];
[maxlen, indx] = max(len, [], 2);

ind0 = [1, 2, 3; 1, 3, 2; 2, 3, 1];

%
fmdl.M_tensor.a = maxlen;

% longest edge vec
lvec = zeros(length(indx),2);
for ii=1:length(indx)
    lvec(ii,:) = fmdl.nodes(fmdl.elems(ii,ind0(indx(ii),2)),:) - fmdl.nodes(fmdl.elems(ii,ind0(indx(ii),1)),:);
end

% abs((y2 - t1)x0 - (x2 - x1)y0 + x2.y1 - y2.x1) / len
ind1 = sub2ind(size(fmdl.elems), (1:size(indx,1)).', ind0(indx, 1));
ind2 = sub2ind(size(fmdl.elems), (1:size(indx,1)).', ind0(indx, 2));
ind3 = sub2ind(size(fmdl.elems), (1:size(indx,1)).', ind0(indx, 3));

fmdl.M_tensor.b = abs((fmdl.nodes(fmdl.elems(ind2),2)- fmdl.nodes(fmdl.elems(ind1),2)).*fmdl.nodes(fmdl.elems(ind3),1) ...
    - (fmdl.nodes(fmdl.elems(ind2),1)- fmdl.nodes(fmdl.elems(ind1),1)).*fmdl.nodes(fmdl.elems(ind3),2) ...
    + fmdl.nodes(fmdl.elems(ind2),2).*fmdl.nodes(fmdl.elems(ind1),1) ...
    - fmdl.nodes(fmdl.elems(ind2),1).*fmdl.nodes(fmdl.elems(ind1),2))./ maxlen;
    
% Angle - a is from x-axis
fmdl.M_tensor.rot = acos(lvec(:,1)./maxlen);




end

