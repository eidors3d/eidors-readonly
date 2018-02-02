function [ fmdl ] = calc_closest_ellipse( fmdl, type )
%CALC_CLOSEST_ELLIPSE 
%   Heuristic closest ellipse, based on longest axis of triangle elements


if nargin==1
    type = 'inellipse';
end

switch type
    
    case 'inellipse'
        % corners
        c1 = fmdl.nodes(fmdl.elems(:,1),:);
        c2 = fmdl.nodes(fmdl.elems(:,2),:);
        c3 = fmdl.nodes(fmdl.elems(:,3),:);
        
        % side lengths
        l1 = sqrt(sum((c2 - c1).^2,2));
        l2 = sqrt(sum((c3 - c1).^2,2));
        l3 = sqrt(sum((c3 - c2).^2,2));
        
        % In Trilinear coords, Steiner inellipse has eqn:
        % (l1 x)^2 + (l2 y)2 + (l3 z)^2 - 2l1 l2 xy - 2 l2 l3 yz - 2 l1 l3 xz = 0
        
        Z = sqrt(l1.^4 + l2.^4 + l3.^4 - l1.^2.*l2.^2 - l2.^2.*l3.^2 - l1.^2.*l3.^2);
        
        fmdl.M_tensor.a = sqrt(l1.^2 + l2.^2 + l3.^2 + 2*Z)/6;
        fmdl.M_tensor.b = sqrt(l1.^2 + l2.^2 + l3.^2 - 2*Z)/6;
        
        % Foci are zeros of derivative of polynomial which has zeros at
        % verticies (Marden's theorem)
        % D (c1(1) + 1i c1(2) - x)(c2(1) + 1i c2(2) - x)(c3(1) + 1i c3(2) -x)=0
        c1 = c1(:,1) + 1i*c1(:,2);
        c2 = c2(:,1) + 1i*c2(:,2);
        c3 = c3(:,1) + 1i*c3(:,2);
        
        a = -3;
        b = 2*(c1 + c2 + c3);
        c = -(c1.*c2 + c1.*c3 + c2.*c3);
        
        % Roots
       st = sqrt(b.^2 - 4*a.*c);
       F1 = (-b + st)./(2*a.*c);
       F1 = [real(F1), imag(F1)];
       
       % Angle between FC and x direction
       FC = F1 - fmdl.elem_centre;
       fmdl.M_tensor.rot = atan(FC(:,2)./FC(:,1));
       
        
        
    case 'heuristic'

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
        fmdl.M_tensor.rot = acos(lvec(:,1)./maxlen) + pi;

end


end

