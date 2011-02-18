function [fwd_model] = linear_reorder(fwd_model,ccw)
%function [fwd_model] = linear_reorder(fwd_model,ccw)
%Function to reorder local nodes (counter)clockwise per element
%Input:  - fwd_model structure
%        - ccw = -1 (default) - counter clockwise OR 1 - clockwise   
%Output: - fwd_model structure (only .elems changes)
%NOTE:Function only for linear triangles, since in this case, identity:
%         No. of nodes/element = No. spatial dimensions + 1

% (C) 2011 Michael Crabb. License: GPL version 2 or version 3
% $Id:$

if (nargin==1) 
    ccw=-1; %Default specify counter-clockwise nodes
end

nodecoords = fwd_model.nodes; %Cache coorindates of nodes [nnodesxnodedim]
elementnodes = fwd_model.elems; %Cache matrix of elements [eletotalxelenode]

eletotal = size(elementnodes,1); %No. of elements
elenode = size(elementnodes,2); %No. of nodes per element

for e=1:eletotal; %Loop over all elements
    %Row vector of global nodes [1xelenode]
    enodes = elementnodes(e,:); 
    %Matrix of nodal positions [elenodexdim] (Linear dimension==elenode-1) 
    nd = nodecoords(enodes,:); 
    
    %Calculate area of triangle/volume defined by the elements nodes
    %In 2D this is area and in 3D this is volume
    area= det([ones(elenode,1),nd]);
    areasign=sign(area); 
    
    %If sign is (pos) neg swap two nodes (last two will suffice..)
    if(areasign == ccw) %Swap last two entries of enodes 
        temp=enodes(elenode-1);
        enodes(elenode-1)=enodes(elenode);
        enodes(elenode) = temp;
        %elementnodes(e,:)=enodes; %Put back into elementnodes matrix
    end
    elementnodes(e,:)=enodes; %Put enodes back into elementnodes matrix
end
fwd_model.elems=elementnodes; %Reassign fwd_model.elems

end %End LinearReorder function
