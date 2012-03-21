function mdl = prefine(mdl)
%PREFINE - Add extra nodes in elements for higher order approximation
%MDL = prefine(MDL)


%PARAMETERS

%Cache the node matrix and element structures and boundary list
nodestruc=mdl.nodes; elemstruc=mdl.elem; boundstruc=mdl.bound;

%Cache no. vertices/boundary edges(faces) and elements
newnodes=size(nodestruc,1); nbounds=size(boundstruc,2); nelems=size(elemstruc,2);

%FUNCTIONS

%Find every elements neighbouring elements (sharing an EDGE)
elemedge = mc_connect_elem_edge_neigh(elemstruc,eletype);

%Find the element to which every boundary belongs
elembound = mc_connect_elem_bound(elemstruc,boundstruc);


%Loop over elems, add new nodes using index vector of element neighbours
% fprintf(1,'Adding extra element nodes\n');
for ii=1:nelems
    %Print out progress
%     if(mod(ii,500)==0); fprintf(1,'%i out of %i elements have been p-refined\n',ii,nelems); end
    
    %2D Element Refinement
    if(strcmp(elemstruc(ii).type,'tri3')) %Do nothing, vertices/nodes equivalent
    elseif(strcmp(elemstruc(ii).type,'tri6'));
        [elemstruc,nodestruc,newnodes] = prefine2dquadelement(elemstruc,nodestruc,newnodes,ii,elemedge); 
    elseif(strcmp(elemstruc(ii).type,'tri10'));
        [elemstruc,nodestruc,newnodes] = prefine2dcubielement(elemstruc,nodestruc,newnodes,ii,elemedge); 
        
    %3D Element Refinement    
    elseif(strcmp(elemstruc(ii).type,'tet4')); %Do nothing, vertices/nodes equivalent
    elseif(strcmp(elemstruc(ii).type,'tet10'));
        [elemstruc,nodestruc,newnodes] = prefine3dquadelement(elemstruc,nodestruc,newnodes,ii,elemedge); 
    end
end


%Go through boundary list and add new nodes (coordinates in node structure
%and equation number in boundstruc (new or old equation number)
% fprintf(1,'Adding extra boundary nodes\n');
for ii=1:nbounds
    %Print out progress
%    if(mod(ii,500)==0); fprintf(1,'%i out of %i boundaries have been p-refined\n',ii,nbounds); end
    %2D Boundary Refinement
    if(strcmp(boundstruc(ii).type,'tri3'))  %Do nothing, vertices/nodes equivalent
    elseif(strcmp(boundstruc(ii).type,'tri6'));
        [boundstruc,elemstruc,nodestruc,newnodes] = prefine2dquadboundary(boundstruc,elemstruc,nodestruc,newnodes,ii,elembound); 
    elseif(strcmp(boundstruc(ii).type,'tri10'));
        [boundstruc,elemstruc,nodestruc,newnodes] = prefine2dcubiboundary(boundstruc,elemstruc,nodestruc,newnodes,ii,elembound); 
        
    %3D Boundary Refinement    
    elseif(strcmp(boundstruc(ii).type,'tet4')); %Do nothing, vertices/nodes equivalent
    elseif(strcmp(boundstruc(ii).type,'tet10'));
        [boundstruc,elemstruc,nodestruc,newnodes] = prefine3dquadboundary(boundstruc,elemstruc,nodestruc,newnodes,ii,elembound); 
    end
end

%Re-assign node, element and boundary structures
mdl.nodes=nodestruc; mdl.elem=elemstruc; mdl.bound=boundstruc;

end





%%%%%%    SPECIFIC FUNCTIONS TO ADD NODES TO BOUNDARY AND ELEMENTS   %%%%%%%





%2D QUADRATIC ELEMENT AND BOUNDARY
function [elemstruc,nodestruc,newnodes] = prefine2dquadelement(elemstruc,nodestruc,newnodes,ii,elemedge)
    
    %Vetex coordinates (by row) of element ii vertices
    thise=nodestruc(elemstruc(ii).nodes,:);
        
    %Find coordinates at midpoints of element ii's edges (new nodes here)
    thise(4,:) = 0.5*(thise(1,:)+thise(2,:)); %12
    thise(5,:) = 0.5*(thise(2,:)+thise(3,:)); %23
    thise(6,:) = 0.5*(thise(1,:)+thise(3,:)); %13
    
    %Find element ii's neighbours (edge neighbours)
    thiseedgeneigh=elemedge(ii).edgeneigh;
    
    %Find number of neighbours
    n_elem_ii_edge_neigh=size(thiseedgeneigh,1);
    
    %Find the first neighbours node numbers
    elemiiedgeneighnode=elemstruc(thiseedgeneigh(1)).nodes;
    
    %Find list of unique node numbers of neighbours
    for pp=1:n_elem_ii_edge_neigh-1
        elemiiedgeneighnodetmp=elemstruc(thiseedgeneigh(pp+1)).nodes;  
        elemiiedgeneighnode=union(elemiiedgeneighnode,elemiiedgeneighnodetmp);
    end

    %Loop over the new nodes
    for j=4:6   
        %List the coordinates of element ii's neighbouring nodes
        nodestrucneigh=nodestruc(elemiiedgeneighnode,:);
        nodestrucdiff=nodestrucneigh;
        
        for i=1:size(nodestrucdiff)
            nodestrucdiff(i,:)=nodestrucneigh(i,:) - thise(j,:);
        end
        %Sum absoluute values across row and test for a 0 index
        sumdifference = sum(abs(nodestrucdiff),2);
        oldnodesloc = find(sumdifference<10^-15);
        
        %Find corresponding global node number
        oldnodes=elemiiedgeneighnode(:,oldnodesloc);
        
        %Test to see if node is already present
        if (size(oldnodes,2)==0)
            %Increase node counter
            newnodes=newnodes+1; 
            %Add coodinates of node to end of node list
            nodestruc(newnodes,:)=thise(j,:);
            %Add the new node number to elemstruc(ii)
            elemstruc(ii).nodes(j)=newnodes; 
            %Add the new node to the neighbour list
            elemiiedgeneighnode=[elemiiedgeneighnode,newnodes];
        else 
            %Add the old node number to elemstruc(ii)
            elemstruc(ii).nodes(j)=oldnodes;
        end
    end
end
function [boundstruc,elemstruc,nodestruc,newnodes] = prefine2dquadboundary(boundstruc,elemstruc,nodestruc,newnodes,ii,elembound)
    
    %Find the vertex coordinates of boundary labelled ii
    thisb=nodestruc(boundstruc(ii).nodes,:);
    
    %Find the coordinate at midpoint of boundary ii
    thisb(3,:) = 0.5*(thisb(1,:)+thisb(2,:));
    
    %Find the node numbers of the element of which boundary belongs
    boundnodenums=elemstruc(elembound(ii)).nodes;
    
    %Find the coordinates of these node numbers
    boundnodecoord=nodestruc(boundnodenums,:);
    
    %Loop over the new nodes
    for j=3       
        %Copy node list and create empty vector of node number
        boundnodecoorddiff=boundnodecoord; 
        
        for i=1:size(boundnodecoorddiff)
            boundnodecoorddiff(i,:)=boundnodecoord(i,:) - thisb(j,:);
        end
        %Sum absoluute values across row, and test if there is a 0 index 
        sumdifference = sum(abs(boundnodecoorddiff),2);
        oldnodesloc = find(sumdifference<10^-15);
        oldnodes=boundnodenums(oldnodesloc);

        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in elem
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1; 
            nodestruc(newnodes,:)=thisb(j,:);
            boundstruc(ii).nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in the elem strucutre
        else
            boundstruc(ii).nodes(j)=oldnodes;
        end
    end
end




%2D CUBIC ELEMENT AND BOUNDARY
function [elemstruc,nodestruc,newnodes] = prefine2dcubielement(elemstruc,nodestruc,newnodes,ii,elemedge)
    %Find the vertex numbers of element
    elemiinodes=elemstruc(ii).nodes;
    
    %Find the vertex coordinate of element labelled ii
    thise=nodestruc(elemiinodes,:);
        
    %Find the coordinates at intervals of 1/3 along edges 
    thise(4,:) = 2.0/3.0*thise(1,:) + 1.0/3.0*thise(2,:);
    thise(5,:) = 1.0/3.0*thise(1,:) + 2.0/3.0*thise(2,:);
    thise(6,:) = 2.0/3.0*thise(2,:) + 1.0/3.0*thise(3,:);
    thise(7,:) = 1.0/3.0*thise(2,:) + 2.0/3.0*thise(3,:);
    thise(8,:) = 2.0/3.0*thise(3,:) + 1.0/3.0*thise(1,:);
    thise(9,:) = 1.0/3.0*thise(3,:) + 2.0/3.0*thise(1,:);
    thise(10,:) = 1.0/3.0*(thise(1,:)+thise(2,:)+thise(3,:));   
    
    %Find element ii's neighbours (edge neighbours)
    thiseedgeneigh=elemedge(ii).edgeneigh;
    
    %Find number of neighbours
    n_elem_ii_edge_neigh=size(thiseedgeneigh,1);
    
    %Find the first neighbours node numbers
    elemiiedgeneighnode=elemstruc(thiseedgeneigh(1)).nodes;
    
    %Find list of unique node numbers of neighbours
    for pp=1:n_elem_ii_edge_neigh-1
        elemiiedgeneighnodetmp=elemstruc(thiseedgeneigh(pp+1)).nodes;  
        elemiiedgeneighnode=union(elemiiedgeneighnode,elemiiedgeneighnodetmp);
    end

    %Loop over the new nodes
    for j=4:10   
        %List the coordinates of element ii's neighbouring nodes
        nodestrucneigh=nodestruc(elemiiedgeneighnode,:);
        nodestrucdiff=nodestrucneigh;
        
        for i=1:size(nodestrucdiff)
            nodestrucdiff(i,:)=nodestrucneigh(i,:) - thise(j,:);
        end
        %Sum absoluute values across row and test for a 0 index
        sumdifference = sum(abs(nodestrucdiff),2);
        oldnodesloc = find(sumdifference<10^-15);
        
        %Find corresponding global node number
        oldnodes=elemiiedgeneighnode(:,oldnodesloc);
        
        %Test to see if node is already present
        if (size(oldnodes,2)==0)
            %Increase node counter
            newnodes=newnodes+1; 
            %Add coodinates of node to end of node list
            nodestruc(newnodes,:)=thise(j,:);
            %Add the new node number to elemstruc(ii)
            elemstruc(ii).nodes(j)=newnodes; 
            %Add the new node to the neighbour list
            elemiiedgeneighnode=[elemiiedgeneighnode,newnodes];
        else 
            %Add the old node number to elemstruc(ii)
            elemstruc(ii).nodes(j)=oldnodes;
        end
    end
end
function [boundstruc,elemstruc,nodestruc,newnodes] = prefine2dcubiboundary(boundstruc,elemstruc,nodestruc,newnodes,ii,elembound)
    
    %Find the vertex coordinates of boundary labelled ii
    thisb=nodestruc(boundstruc(ii).nodes,:);
    
    %Find the coordinate at intervals of 1/3, 2/3 along boundary ii
    thisb(3,:) = 2.0/3.0*thisb(1,:)+1.0/3.0*thisb(2,:);
    thisb(4,:) = 1.0/3.0*thisb(1,:)+2.0/3.0*thisb(2,:);
    
    %Find the node numbers of the element of which boundary belongs
    boundnodenums=elemstruc(elembound(ii)).nodes;
    
    %Find the coordinates of these node numbers
    boundnodecoord=nodestruc(boundnodenums,:);
    
    %Loop over the new nodes
    for j=3:4       
        %Copy node list and create empty vector of node number
        boundnodecoorddiff=boundnodecoord; 
        
        for i=1:size(boundnodecoorddiff)
            boundnodecoorddiff(i,:)=boundnodecoord(i,:) - thisb(j,:);
        end
        %Sum absoluute values across row, and test if there is a 0 index 
        sumdifference = sum(abs(boundnodecoorddiff),2);
        oldnodesloc = find(sumdifference<10^-15);
        oldnodes=boundnodenums(oldnodesloc);

        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in elem
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1; 
            nodestruc(newnodes,:)=thisb(j,:);
            boundstruc(ii).nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in the elem strucutre
        else
            boundstruc(ii).nodes(j)=oldnodes;
        end
    end
end



%3D QUADRATIC TETRAHEDRA ELEMENT AND BOUNDARY
function [elemstruc,nodestruc,newnodes] = prefine3dquadelement(elemstruc,nodestruc,newnodes,ii,elemedge)

    %Find the vertex coordinates corresponding to element ii
    thise=nodestruc(elemstruc(ii).nodes,:);
        
    %Calculate new coordinates (in correct order to match with the 
    %positioning of shape functions) at midpoints
    thise(5,:) = 0.5*(thise(1,:)+thise(2,:)); %12
    thise(6,:) = 0.5*(thise(1,:)+thise(3,:)); %13
    thise(7,:) = 0.5*(thise(1,:)+thise(4,:)); %14
    thise(8,:) = 0.5*(thise(2,:)+thise(3,:)); %23
    thise(9,:) = 0.5*(thise(3,:)+thise(4,:)); %34
    thise(10,:) = 0.5*(thise(2,:)+thise(4,:)); %24
        
    %Find element ii's neighbours (edge neighbours)
    thiseedgeneigh=elemedge(ii).edgeneigh;
    
    %Find number of neighbours
    n_elem_ii_edge_neigh=size(thiseedgeneigh,1);
    
    %Find the first neighbours node numbers
    elemiiedgeneighnode=elemstruc(thiseedgeneigh(1)).nodes;
    
    %Find list of unique node numbers of neighbours
    for pp=1:n_elem_ii_edge_neigh-1
        elemiiedgeneighnodetmp=elemstruc(thiseedgeneigh(pp+1)).nodes;  
        elemiiedgeneighnode=union(elemiiedgeneighnode,elemiiedgeneighnodetmp);
    end

    %Loop over the new nodes
    for j=4:10   
        %List the coordinates of element ii's neighbouring nodes
        nodestrucneigh=nodestruc(elemiiedgeneighnode,:);
        nodestrucdiff=nodestrucneigh;
        
        for i=1:size(nodestrucdiff)
            nodestrucdiff(i,:)=nodestrucneigh(i,:) - thise(j,:);
        end
        %Sum absoluute values across row and test for a 0 index
        sumdifference = sum(abs(nodestrucdiff),2);
        oldnodesloc = find(sumdifference<10^-15);
        
        %Find corresponding global node number
        oldnodes=elemiiedgeneighnode(:,oldnodesloc);
        
        %Test to see if node is already present
        if (size(oldnodes,2)==0)
            %Increase node counter
            newnodes=newnodes+1; 
            %Add coodinates of node to end of node list
            nodestruc(newnodes,:)=thise(j,:);
            %Add the new node number to elemstruc(ii)
            elemstruc(ii).nodes(j)=newnodes; 
            %Add the new node to the neighbour list
            elemiiedgeneighnode=[elemiiedgeneighnode,newnodes];
        else 
            %Add the old node number to elemstruc(ii)
            elemstruc(ii).nodes(j)=oldnodes;
        end
    end
end
function [boundstruc,elemstruc,nodestruc,newnodes] = prefine3dquadboundary(boundstruc,elemstruc,nodestruc,newnodes,ii,elembound)
    
    %Find the vertex coordinate of boundary labelled ii
    thisb=nodestruc(boundstruc(ii).nodes,:);
        
    %Find the coordinates at midpoints of element ii's edges
    thisb(4,:) = 0.5*(thisb(1,:)+thisb(2,:)); %12
    thisb(5,:) = 0.5*(thisb(1,:)+thisb(3,:)); %13
    thisb(6,:) = 0.5*(thisb(2,:)+thisb(3,:)); %23
    
    %Find the node numbers of the element of which boundary belongs
    boundnodenums=elemstruc(elembound(ii)).nodes;
    
    %Find the coordinates of these node numbers
    boundnodecoord=nodestruc(boundnodenums,:);
    
    %Loop over the new nodes
    for j=3:4       
        %Copy node list and create empty vector of node number
        boundnodecoorddiff=boundnodecoord; 
        
        for i=1:size(boundnodecoorddiff)
            boundnodecoorddiff(i,:)=boundnodecoord(i,:) - thisb(j,:);
        end
        %Sum absoluute values across row, and test if there is a 0 index 
        sumdifference = sum(abs(boundnodecoorddiff),2);
        oldnodesloc = find(sumdifference<10^-15);
        oldnodes=boundnodenums(oldnodesloc);

        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in elem
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1; 
            nodestruc(newnodes,:)=thisb(j,:);
            boundstruc(ii).nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in the elem strucutre
        else
            boundstruc(ii).nodes(j)=oldnodes;
        end
    end
end


%OLD CODE (WITHOUT INDEXING ARRAYS OF NEIGHBOURING EDGES)

%2D CUBIC ELEMENT AND BOUNDARY (OLD CODE)
%{
function [elemstrucii,nodestruc,newnodes] = prefine2dcubielement(elemstrucii,nodestruc,newnodes)
    %Find the vertex numbers of element
    elemiinodes=elemstrucii.nodes;
    
    %Find the vertex coordinate of element labelled ii
    thise=nodestruc(elemiinodes,:);
        
    %Find the coordinates at intervals of 1/3 along edges 
    thise(4,:) = 2.0/3.0*thise(1,:) + 1.0/3.0*thise(2,:);
    thise(5,:) = 1.0/3.0*thise(1,:) + 2.0/3.0*thise(2,:);
    thise(6,:) = 2.0/3.0*thise(2,:) + 1.0/3.0*thise(3,:);
    thise(7,:) = 1.0/3.0*thise(2,:) + 2.0/3.0*thise(3,:);
    thise(8,:) = 2.0/3.0*thise(3,:) + 1.0/3.0*thise(1,:);
    thise(9,:) = 1.0/3.0*thise(3,:) + 2.0/3.0*thise(1,:);
    thise(10,:) = 1.0/3.0*(thise(1,:)+thise(2,:)+thise(3,:));   
    
    
    
    %Loop over new nodes, find pointwise difference of coordinates with all 
    %entries in the node matrix. Two possibilities:
    %1. Pointwise difference is 0 for 1 index only (I think there can only
    %be one such index....), so assign the known index in correct local 
    %position in mdl.elem
    %2. Pointwise difference is NOT 0, assign a new node number, add the
    %coordinates in node matrix and put the node number in correct local
    %position in mdl.elem
    
    for j=4:10   
        %Copy node list and create empty vector of node number
        nodestrucdiff=nodestruc; 
        
        for i=1:size(nodestruc)
            nodestrucdiff(i,:)=nodestruc(i,:) - thise(j,:);
        end
        %Sum absoluute values across row, and test if there is a 0 index 
        sumdifference = sum(abs(nodestrucdiff),2);
        oldnodes = find(sumdifference<10^-15); %TOLERANCE NEEDED
        
        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in elem
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1; 
            nodestruc(newnodes,:)=thise(j,:);
            elemstrucii.nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in the elem strucutre
        else
            elemstrucii.nodes(j)=oldnodes;
        end
    end
end

function [boundstrucii,nodestruc,newnodes] = prefine2dcubiboundary(boundstrucii,nodestruc,newnodes)
    
    %Find the vertex coordinates of boundary labelled ii
    thisb=nodestruc(boundstrucii.nodes,:);
    
    %Find the coordinate at intervals of 1/3, 2/3 along boundary ii
    thisb(3,:) = 2.0/3.0*thisb(1,:)+1.0/3.0*thisb(2,:);
    thisb(4,:) = 1.0/3.0*thisb(1,:)+2.0/3.0*thisb(2,:);
    
    %1. Unique node number for additional node required on boundary
    %Loop over new nodes, find pointwise difference of coordinates with all 
    %entries in the node matrix. Two possibilities:
    %1. Pointwise difference is 0 for 1 index only (I think there can only
    %be one such index....), so assign the known index in correct local 
    %position in mdl.bound
    %2. Pointwise difference is NOT 0, assign a new node number, add the
    %coordinates in node matrix and put the node number in correct local
    %position in mdl.bound
    
    for j=3:4       
        %Copy node list and create empty vector of node number
        nodestrucdiff=nodestruc; 
        
        for i=1:size(nodestruc)
            nodestrucdiff(i,:)=nodestruc(i,:) - thisb(j,:);
        end
        %Sum absoluute values across row, and test if there is a 0 index 
        sumdifference = sum(abs(nodestrucdiff),2);
        oldnodes = find(sumdifference<10^-15);

        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in elem
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1; 
            nodestruc(newnodes,:)=thisb(j,:);
            boundstrucii.nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in the elem strucutre
        else
            boundstrucii.nodes(j)=oldnodes;
        end
    end
end
%}

%3D QUADRATIC TETRAHEDRA ELEMENT AND BOUNDARY (OLD CODE)
%{
function [boundstrucii,nodestruc,newnodes] = prefine3dquadboundary(boundstrucii,nodestruc,newnodes)
    
    %Find the vertex coordinate of boundary labelled ii
    thisb=nodestruc(boundstrucii.nodes,:);
        
    %Find the coordinates at midpoints of element ii's edges
    thisb(4,:) = 0.5*(thisb(1,:)+thisb(2,:)); %12
    thisb(5,:) = 0.5*(thisb(1,:)+thisb(3,:)); %13
    thisb(6,:) = 0.5*(thisb(2,:)+thisb(3,:)); %23
    
    %Loop over new nodes, find pointwise difference of coordinates with all 
    %entries in the node matrix. Two possibilities:
    %1. Pointwise difference is 0 for 1 index only (I think there can only
    %be one such index....), so assign the known index in correct local 
    %position in mdl.elem
    %2. Pointwise difference is NOT 0, assign a new node number, add the
    %coordinates in node matrix and put the node number in correct local
    %position in mdl.elem
    for j=4:6       
        %Copy node list and create empty vector of node number
        nodestrucdiff=nodestruc;
        
        %Faster Code
        if 1 %Old code
            for i=1:size(nodestruc)
               nodestrucdiff(i,:)=nodestruc(i,:) - thisb(j,:);
            end
            %Sum absoluute values across row, and test if there is a 0 index 
            sumdifference = sum(abs(nodestrucdiff),2);
            oldnodes = find(sumdifference<10^-15); %TOLERANCE NEEDED
        end
        
        %Alternative slower code
        if 0
            oldnodes=[]; %Create an empty oldnodes vector
            %Find pointwise difference of the new node coordinate with old
            for i=1:size(nodestruc) %New code
                nodestrucdiff(i,:)=nodestruc(i,:) - thisb(j,:);
                sumdiff = sum(abs(nodestrucdiff(i)),2);
                if(sumdiff < 10^-15 )
                    oldnodes = i;
                break; %Exit for loop (only 1 possible index)
                end
            end
        end
        
        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in bound
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1; 
            nodestruc(newnodes,:)=thisb(j,:);
            boundstrucii.nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in bound strucutre
        else
            boundstrucii.nodes(j)=oldnodes;
        end
    end
end
function [elemstrucii,nodestruc,newnodes] = prefine3dquadelement(elemstrucii,nodestruc,newnodes)

    %Find the vertex coordinates corresponding to element ii
    thise=nodestruc(elemstrucii.nodes,:);
        
    %Calculate new coordinates (in correct order to match with the 
    %positioning of shape functions) at midpoints
    thise(5,:) = 0.5*(thise(1,:)+thise(2,:)); %12
    thise(6,:) = 0.5*(thise(1,:)+thise(3,:)); %13
    thise(7,:) = 0.5*(thise(1,:)+thise(4,:)); %14
    thise(8,:) = 0.5*(thise(2,:)+thise(3,:)); %23
    thise(9,:) = 0.5*(thise(3,:)+thise(4,:)); %34
    thise(10,:) = 0.5*(thise(2,:)+thise(4,:)); %24
        
    %Loop over new nodes, find pointwise difference of coordinates with all 
    %entries in the node matrix. Two possibilities:
    %1. Pointwise difference is 0 for 1 index only (I think there can only
    %be one such index....), so assign the known index in correct local 
    %position in mdl.elem
    %2. Pointwise difference is NOT 0, assign a new node number, add the
    %coordinates in node matrix and put the node number in correct local
    %position in mdl.elem
    for j=5:10      
        %Copy node list and create empty vector of node number
        nodestrucdiff=nodestruc;  
        
        %Faster Code
        if 1 
            for i=1:size(nodestruc)
                nodestrucdiff(i,:)=nodestruc(i,:) - thise(j,:);
            end
            %Sum absoluute values across row, and test if there is a 0 index 
            sumdifference = sum(abs(nodestrucdiff),2);
            oldnodes = find(sumdifference<10^-15); %TOLERANCE NEEDED
        end
        
        %Alternative slower code
        if 0 
            oldnodes=[]; %Create empty old nodes vector 
            %Find pointwise difference of the new node coordinate with old
            for i=1:size(nodestruc) %New code
                nodestrucdiff(i,:)=nodestruc(i,:) - thise(j,:);
                sumdiff = sum(abs(nodestrucdiff(i)),2);
                if(sumdiff < 10^-15 )
                    oldnodes = i;
                    break; %Exit for loop (only 1 possible index)
                end
            end
        end
        
        %If no index, then this is a new node, assign new equation number,
        %put coords in node structure and new equation number in elem
        if (size(oldnodes,1)==0)
            newnodes=newnodes+1;
            nodestruc(newnodes,:)=thise(j,:);
            elemstrucii.nodes(j)=newnodes; 
            
        %If index is here (only 1 such index), put the old equation number
        %in the correct local position in the elem strucutre
        else
            elemstrucii.nodes(j)=oldnodes;
        end
    end
end
%}
end
