function [mdl] = linearelemreorder(mdl,ccw)
%function [mdl] = linearelemreorder(mdl,ccw)
%Function to reorder local vertices (counter)clockwise in element 
%Input:  - mdl structure
%        - ccw = -1 (default) - counter clockwise OR 1 - clockwise   
%Output: - mdl structure (only mdl.elem(i) changes)

if (nargin==1) 
    ccw=-1;  %Default specify : Counter-clockwise relative to first node in 2D
	     %                : Clockwise looking INTO relative to first node in 3D
end

%Cache node matrix and calculate number of elements
nodecoords = mdl.nodes; nelems=size(mdl.elem,2); 

for e=1:nelems;
    %Row vector of the node numbers and the no. of vertices
    enodes = mdl.elem(e).nodes; elenode = size(enodes,2);
    
    %Matrix of nodal positions [elenodexdim] (Linear dimension==elenode-1) 
    nd = nodecoords(enodes,:);
    
    %Calculate area(2D)/volume(3D) defined by the vertices
    area= det([ones(elenode,1),nd]); areasign=sign(area);
    
    %If sign is (pos) neg swap two nodes (last two will suffice..)
    if(areasign == ccw) %Swap last two entries of enodes 
        temp=enodes(elenode-1);
        enodes(elenode-1)=enodes(elenode);
        enodes(elenode) = temp;
    end
    mdl.elem(e).nodes=enodes; %Put enodes back into elementnodes matrix
end

end
