function [bound,elem,nodes]=fem_1st_to_higher_order(fwd_model)
% FEM_MODIFY:  Modify the FEM for high order FEM called as
%    [bound,elem] = mc_fem_modify( fwd_model )
% It will call :fwd_model.modify=@mc_fem_modify
% where we need to make sure fwd_model.approx_type='tri3', 'tet4' etc. 
%
% fwd_model : is a fwd_model structure
% bound : is the boundary nodes numbers (from boundary (bound(ii).nodes))
% elem : is the element nodes numbers (from elems elem(ii).nodes))
%
%M Crabb - 29.06.2012

bound = eidors_obj('get-cache', {fwd_model}, 'bound');
elem = eidors_obj('get-cache', {fwd_model}, 'elem');
nodes=eidors_obj('get-cache',{fwd_model},'nodes');
if ~isempty(bound)
    if ~isempty(elem)
        if ~isempty(nodes)
            eidors_msg('bound: using cached value', 3);        
            eidors_msg('elem: using cached value',3);
            eidors_msg('nodes: using cached value',3);
        return
        end
    end
end

[bound,elem,nodes]= mc_fem_modify(fwd_model);

eidors_obj('set-cache', {fwd_model}, 'bound', bound);
eidors_obj('set-cache', {fwd_model}, 'elem',elem);
eidors_obj('set-cache', {fwd_model}, 'nodes',nodes);
eidors_msg('fem_modify: setting cached value', 3);
end

function [bound,elem,nodes]=mc_fem_modify(fwd_model)
%Function  that takes a fwd_model structure (with at least fwd_model.boundary,
%fwd_model.elems, fwd_model.nodes), reorders these to form anti-clockwise sets, and
%then adds extra nodes for p-refinement.
%
%Function exits with fwd_model.elem(i).nodes; fwd_model.bound(i).nodes, which have the
%extra p-refined nodes correctly positioned, and still has the fwd_model.boundary
%and fwd_model.elems structures

%Change to fwd_model.elem(i) - do from last to first
for i=size(fwd_model.elems,1):-1:1
    fwd_model.elem(i).nodes=fwd_model.elems(i,:);
end

%Change to fwd_model.boundary(i) - do from last to first
for i=size(fwd_model.boundary,1):-1:1
    fwd_model.bound(i).nodes=fwd_model.boundary(i,:);
end

%Prefine - boundary first, then remaining elements
fwd_model=mc_p_refine(fwd_model);

%FIXME Clear this up, by changing the output of two functions above, 
%so that the output is just bound and elem...
bound=fwd_model.bound; elem=fwd_model.elem; nodes=fwd_model.nodes;
end


function mdl = mc_p_refine(mdl)
%PROPER WAY TO DO THIS IN 2D
%1. CALCULATE UNIQUE EDGES THROUGH THE TOPOLOGY ARRAY (ELEMS)
%2. CALCULATE DUAL GRAPH THROUGH THE EDGES
%3. ADD NODES BY TRAVERSING THROUGH THE DUAL GRAPH 
%4. INSTEAD OF TESTING ALL NODES FIND THE NEIGBOURING ELEMENTS BY EDGES BEFORE P-REFINE
%PREFINE - Add extra nodes in elements for higher order approximation
%MDL = mc_p_refine(MDL)

%Find the node/elem/bound structure
nodestruc=mdl.nodes; boundstruc = mdl.bound; elemstruc = mdl.elem;

%Get the no. nodes/bounds/elems and dimension
nodedim=size(nodestruc,2); newnodes=size(nodestruc,1); nbounds=size(boundstruc,2); nelems=size(elemstruc,2);

%Local ordering orientation parameter 
ccw=-1;

%Reorder element nodes (locally) anti-clockwise relative to first node
elemstruc = linear_elem_reorder(elemstruc,nodestruc,ccw);

%Find vector (n_bounds) of which element is ith boundary a member of
elembound = mc_connect_elem_bound(elemstruc,boundstruc,nodedim);

%Reorder boundary nodes (locally) anti-clockwise relative to interior node
boundstruc = linear_bound_reorder(boundstruc,nodestruc,elemstruc,elembound,ccw);

%Find every elements neighbouring elements (sharing an EDGE)
elemedge = mc_connect_elem_edge_neigh(elemstruc,mdl.approx_type,nodedim);

%Loop over elems, add new nodes using index vector of element neighbours
if(strcmp(mdl.approx_type,'tri3') || strcmp(mdl.approx_type,'tet4'))
    %Do nothing, need no refinement
elseif(strcmp(mdl.approx_type,'tri6'))    
    for ii=1:nelems
        [elemstruc,nodestruc,newnodes] = prefine2dquadelement(elemstruc,nodestruc,newnodes,ii,elemedge); 
    end
    for jj=1:nbounds
        [boundstruc,elemstruc,nodestruc,newnodes] = prefine2dquadboundary(boundstruc,elemstruc,nodestruc,newnodes,jj,elembound); 
    end        
elseif(strcmp(mdl.approx_type,'tri10'))
    for ii=1:nelems
        [elemstruc,nodestruc,newnodes] = prefine2dcubielement(elemstruc,nodestruc,newnodes,ii,elemedge); 
    end
    for jj=1:nbounds
        [boundstruc,elemstruc,nodestruc,newnodes] = prefine2dcubiboundary(boundstruc,elemstruc,nodestruc,newnodes,jj,elembound); 
    end
elseif(strcmp(mdl.approx_type,'tet10'));
    for ii=1:nelems
        [elemstruc,nodestruc,newnodes] = prefine3dquadelement(elemstruc,nodestruc,newnodes,ii,elemedge); 
    end
    for jj=1:nbounds
         [boundstruc,elemstruc,nodestruc,newnodes] = prefine3dquadboundary(boundstruc,elemstruc,nodestruc,newnodes,jj,elembound);        
    end
else
    error('mc_fem_modify: Element type ("%s") not recognised',mdl.approx_type);
end

%Re-assign node, element and boundary structures
mdl.nodes=nodestruc; mdl.elem=elemstruc; mdl.bound=boundstruc;

end



%SPECIFIC FUNCTIONS HERE FOR THE FUNCTION


%1.GEOMETRY FUNCTIONS

%(a) FUNCTION TO FIND NEIGHBOURING ELEMENTS (BY EDGES) OF EACH ELEMENT
function [elemedge]=mc_connect_elem_edge_neigh(elemstruc,eletype,nodedim)

    %Print out progress
%    tic; fprintf(1,'Connect each element with its neighbouring elements (edges) ');
    if (strcmp(eletype,'tri3') || strcmp(eletype,'tet4'))
        %Do nothing - prefine.m redundant for lienar approx don't need neighbours 
        elemedge=[];
    else
        %Find no. of elements and revert to old matrix structure
        nelems=size(elemstruc,2);
        for j=1:nelems
            elemstrucold(j,:)=elemstruc(j).nodes;
        end
        %Loop over elements
        for ii=1:nelems
            %Vector of vertex numbers of this boundary
            elemiinodes=elemstruc(ii).nodes;
    
            %TODO _ SPEED THIS UP!!!!!
            if(nodedim==2) %2D problemm
                %Find elements to which each node of elem ii belongs
                [rownode1,blah]=find(elemstrucold==elemiinodes(1));
                [rownode2,blah]=find(elemstrucold==elemiinodes(2));
                [rownode3,blah]=find(elemstrucold==elemiinodes(3));
    
                %Intersection vectors above gives common elements along edge
                r1r2=intersect(rownode1,rownode2); %edge 1
                r1r3=intersect(rownode1,rownode3); %edge 2
                r2r3=intersect(rownode2,rownode3); %edge 3     
        
                %Now we find the union of the vectors
                union1=union(r1r2,r1r3);
                union2=union(union1,r2r3);
            
                %Now eliminate the current element from this 
                elemiiedgeneighbour=setdiff(union2,ii);
        
                %Store this in mdl.elem(ii)
                elemedge(ii).edgeneigh=elemiiedgeneighbour;                
            elseif(nodedim==3) %3D problem
                %Find elements to which each node of elem ii belongs
                [rownode1,blah]=find(elemstrucold==elemiinodes(1));
                [rownode2,blah]=find(elemstrucold==elemiinodes(2));
                [rownode3,blah]=find(elemstrucold==elemiinodes(3)); 
                [rownode4,blah]=find(elemstrucold==elemiinodes(4));
        
                %Intersection vectors above gives common elements along edge
                r1r2=intersect(rownode1,rownode2);
                r1r3=intersect(rownode1,rownode3);
                r1r4=intersect(rownode1,rownode4);
                r2r3=intersect(rownode2,rownode3);
                r2r4=intersect(rownode2,rownode4);
                r3r4=intersect(rownode3,rownode4);
        
                %Now we find the union of the vectors
                union1=union(r1r2,r1r3);
                union2=union(r1r4,r2r3);
                union3=union(r2r4,r3r4);
                union4=union(union1,union2);
                union5=union(union3,union4);
            
                %Now eliminate the current element from this 
                elemiiedgeneighbour=setdiff(union5,ii);
        
                %Store this in mdl.elem(ii)
                elemedge(ii).edgeneigh=elemiiedgeneighbour;
            end 
        end     
    end
% toc;
end

%(b) FUNCTION TO FIND EACH BOUNDARY (UNIQUE) MEMBER ELEMENT
function [bound_elem]=mc_connect_elem_bound(elemstruc,boundstruc,nodedim)
    %Find no. of elements and revert to old matrix structure
    nelems=size(elemstruc,2);
    for j=1:nelems
        elemstrucold(j,:)=elemstruc(j).nodes;
    end
    %Find no. of boundary
    nbounds=size(boundstruc,2);
    for ii=1:nbounds
        %Vector of vertes numbers of this boundary
        boundaryiinodes=boundstruc(ii).nodes;
        %OLD CODE
        if 1
            if(nodedim==2) %2D problem
                %Find row(s) of elems for which each boundary vertex belongs
                [rownode1,blah]=find(elemstrucold==boundaryiinodes(1));
                [rownode2,blah]=find(elemstrucold==boundaryiinodes(2));
        
                %Intersection of rownode1 and rownode2 is the unique element
                boundiielem=intersect(rownode1,rownode2);
            elseif(nodedim==3) %3D problem
                %Find row(s) of elems for which each node belongs
                [rownode1,blah]=find(elemstrucold==boundaryiinodes(1));
                [rownode2,blah]=find(elemstrucold==boundaryiinodes(2));
                [rownode3,blah]=find(elemstrucold==boundaryiinodes(3));

                %Intersection of rownode1 and rownode2 gives a choice
                rownode1node2=intersect(rownode1,rownode2);

                %Intersection of rownode3 with vector above is unique element
                boundiielem=intersect(rownode3,rownode1node2);  
            end
        end
        %NEW CODE (NOT WORKING YET)
        if 0
        if(nodedim==2) %2D problem
            %Find row(s) of elems of which first boundary vertex belongs
            [rownode1,blah]=find(elemstrucold==boundaryiinodes(1));
            %Find row(s) of elems for the boundaries above
            [rownode2,blah]=find(elemstrucold(rownode1,:)==boundaryiinodes(2));
        
            %This is now our unique? element
            boundiielem=rownode1(rownode2);
            
        elseif(nodedim==3) %3D problem
            %Find row(s) of elems for which each node belongs
            [rownode1,blah]=find(elemstrucold==boundaryiinodes(1));
            
            %Find row(s) of elems for the boundaries above
            [rownode2,blah]=find(elemstrucold(rownode1,:)==boundaryiinodes(2));
            rownode2=rownode1(rownode2);
            
             %Find row(s) of elems for the boundaries above
            [rownode3,blah]=find(elemstrucold(rownode2,:)==boundaryiinodes(3));

            %Intersection of rownode3 with vector above is unique element
            boundiielem=rownode2(rownode3);  
        end
        end
        %Store this unique number in a vector
        %FIXME!!! 3D Inclusion models seem to include the boundaries of the
        %inclusions as boundaries. In this case there may be multiple
        %boundaries. Just choose first one for time being, and if this
        %is ACTUAL boundary, then should be unique.....
        bound_elem(ii)=boundiielem(1);
    end    
end



%2. LOCAL RENUMBERING OF NODES (TO GET CLCOKWISE CONVENTION
%STEP 1 : Reorder nodes locally within the element, and boundary to get a
%consistent oriented element
%Elem Default : Arrange nodes locally in mdl.elem so that:
%2D: Counter-clockwise looking down onto element
%3D: Counter-clockwise looking into element from first node
%    (clockwise looking down from ouside element opposite first node)
%
%Boundary Default: Arrange nodes locally in mdl.bound so that
%2D: Counter-clockwise looking down relative to interior node
%3D: Clockwise relative looking into element from interior node
%    (counter-clockwise looking down from outside element opposite interior node)
%
%These conventions are chosen to match the reference element geometry

%(a)REORDERIJNG ELEMENT NUMBERS
function [elemstruc] = linear_elem_reorder(elemstruc,nodestruc,ccw)
    nelems=size(elemstruc,2);
    for e=1:nelems;
        %Row vector of the node numbers and the no. of vertices
        enodes = elemstruc(e).nodes; elenode = size(enodes,2);
    
        %Matrix of nodal positions [elenodexdim] (Linear dimension==elenode-1) 
        nd = nodestruc(enodes,:);
    
       %Calculate area(2D)/volume(3D) defined by the vertices
        area= det([ones(elenode,1),nd]); areasign=sign(area);
    
        %If sign is (pos) neg swap two nodes (last two will suffice..)
        if(areasign == ccw) %Swap last two entries of enodes 
            temp=enodes(elenode-1);
            enodes(elenode-1)=enodes(elenode);
            enodes(elenode) = temp;
        end
        elemstruc(e).nodes=enodes; %Put enodes back into elementnodes matrix
    end
end

%REORDERING BOUNDARY NUMBERS
function [boundstruc]=linear_bound_reorder(boundstruc,nodestruc,elemstruc,bound_elem,ccw)
    nbounds=size(boundstruc,2);
    for bb=1:nbounds;    
        %Find the global node numbers associated with boundary
        bnodes=boundstruc(bb).nodes; 
    
        %Find the node numbers of the element on which boundary is on
        beleind=bound_elem(bb); belenodes=elemstruc(beleind).nodes;
  
        %Find unique internal node and order element node so internal is first 
        intnode=setdiff(belenodes,bnodes); elemboundnode=[intnode,bnodes];
    
        %List by row coordinates of the element
        nodepos=nodestruc(elemboundnode,:); nvertices=size(belenodes,2);     

        %Calculate area(2D)/volume(3D) defined by the vertices
        area= det([ones(nvertices,1),nodepos]); areasign=sign(area);
    
        %If positive area, swap the last two nodes
        if(areasign == -ccw) %Swap last two entries of enodes 
            temp=elemboundnode(end-1);
            elemboundnode(end-1)=elemboundnode(end);
            elemboundnode(end) = temp;
        end
    
        %Put the node numbers back into mdl.bound(bb) (excluding internal node)
        boundstruc(bb).nodes=elemboundnode(2:end);
    end          
end



%3. SPECIFIC FUNCTIONS TO ADD NODES TO BOUNDARY AND ELEMENTS



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
    for j=4:6       
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
