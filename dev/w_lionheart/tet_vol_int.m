function volint=tet_vol_int(v1,v2)
% Find vertices of intersection volume of two tetrahedra

% List of choices of 2 from 1 and one from the other
choices =[ 1,2,1;1,2,2;1,2,3;1,2,4;
           1,3,1;1,3,2;1,3,3;1,3,4;
           1,4,1;1,4,2;1,4,3;1,4,4;                      
           2,3,1;2,3,2;2,3,3;2,3,4;
           2,4,1;2,4,2;2,4,3;2,4,4;
           3,4,1;3,4,2;3,4,3;3,4,4];
           
           
           
    
epsilon=1e-10;
[A1,b1]=tet_to_inequal(v1);
[A2,b2]=tet_to_inequal(v2);
i12 = (A1*v2' -b1*ones(1,4))< epsilon;
i21 = (A2*v1' -b2*ones(1,4))< epsilon;
if all(all(i12))
    volint=tet_vol(v2);
elseif all(all(i21))
    volint=tet_vol(v1);
elseif all(~all(i12)) & all(~all(i21));
    % disjoint
    volint=0;
else
% some intersection, 
vs=[];
% add the vertices that are already in both
vs=[vs;v1(find(all(i21)),:)];
vs=[vs;v2(find(all(i12)),:)];

%try all two faces from one intersected with one face from the
%other
 for i = 1:24
    A=[A1(choices(i,1:2),:);A2(choices(i,3),:)];
    b=[b1(choices(i,1:2));b2(choices(i,3))];
    
    if abs(det(A))>1e-5
     vs=[vs;(A\b)'];
    end
      
    A=[A2(choices(i,1:2),:);A1(choices(i,3),:)];
    b=[b2(choices(i,1:2));b1(choices(i,3))];
 
    if abs(det(A))>1e-5
     vs=[vs;(A\b)'];
    end
    
 end
 if isempty(vs)
    volint=0;
 else
   thosein = find( all(A1 * vs'- b1*ones(1,size(vs,1),1) <epsilon) &  all(A2 * vs'- b2*ones(1,size(vs,1)) <epsilon,1));
   vs=vs(thosein,:);
   vs=unique(vs,'rows');
   if size(vs,1)<4
       volint=0;
   else   
       [K,volint]=convhulln(vs);
   end
 end    
end
    
end

function vol=tet_vol(v)
%finds the volume of one tetrahedron
edges= v(2:end,:)-ones(3,1)*v(1,:);
vol= abs(det(edges))/6;
end   
    
