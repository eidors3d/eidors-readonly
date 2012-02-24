function mdl = mc_elem_bound_reorder(mdl)
%Elem Default : Arrange nodes locally in mdl.elem so that:
%2D: Counter-clockwise looking down onto element
%3D: Clockwise looking into element from first node
%
%Boundary Default: Arrange nodes locally in mdl.bound so that
%2D: Counter-clockwise looking down relative to interior node
%3D: Clockwise relative looking into element from interior node
%
%These conventions are chosen to match the reference element geometry


%PARAMETERS

%Find node/elem/bound structure and problem dimension
nodestruc=mdl.nodes; nodedim=size(nodestruc,2); boundstruc = mdl.bound; elemstruc = mdl.elem;

%This specified correct arrangement as above
ccw=-1;


%FUNCTIONS

%Reorder element nodes (locally) anti-clockwise relative to first node
elemstruc = linear_elem_reorder(elemstruc,nodestruc);

%Find vector (n_bounds) of which element is ith boundary a member of
bound_elem = mc_connect_elem_bound(elemstruc,boundstruc,nodedim);

%Reorder boundary nodes (locally) anti-clockwise relative to interior node
boundstruc = linear_bound_reorder(boundstruc,nodestruc,bound_elem);



%Store the elem and bound structure back in the model
mdl.elem = elemstruc; mdl.bound = boundstruc;



function [elemstruc] = linear_elem_reorder(elemstruc,nodestruc)
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


function [boundstruc]=linear_bound_reorder(boundstruc,nodestruc,bound_elem)
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
    
        %If negative area, swap the 
        if(areasign == ccw) %Swap last two entries of enodes 
            temp=elemboundnode(end-1);
            elemboundnode(end-1)=elemboundnode(end);
            elemboundnode(end) = temp;
        end
    
        %Put the node numbers back into mdl.bound(bb) (excluding internal node)
        boundstruc(bb).nodes=elemboundnode(2:end);
    end          
end

%Put elemstruc and boundstruc back into mdl
mdl.elem=elemstruc; mdl.bound=boundstruc;

end