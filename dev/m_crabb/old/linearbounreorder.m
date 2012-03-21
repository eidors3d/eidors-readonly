function [mdl] = linearbounreorder(mdl,ccw)
%function [mdl] = linearbounreorder(mdl,ccw)
%Function to reorder local vertices (counter)clockwise in boundary relative 
%to the internal node of boundary's element
%
%Input:  - mdl structure
%        - ccw = -1 (default) - counter clockwise OR 1 - clockwise   
%Output: - mdl structure (only mdl.bound(i) changes)

if (nargin==1) 
    ccw=-1; %Default specify : Counter-clockwise relative to interior node in 2D
	    %                : Clockwise looking INTO relative to interior node in 3D
end

%Find node structure and the problem dimension 
nodestruc = mdl.nodes;

%Find boundary sturcture and the no. of boundaries
boundstruc = mdl.bound; nbounds=size(boundstruc,2);

%Find elem structure and the no. of elements
elemstruc=mdl.elem; 

%Loop over the boundaries
for bb=1:nbounds;    
    %Find the global node numbers associated with boundary
    bnodes=boundstruc(bb).nodes; 
    
    %Find the node numbers of the element on which boundary is on
    beleind=boundstruc(bb).elem; belenodes=elemstruc(beleind).nodes;
  
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
    mdl.bound(bb).nodes=elemboundnode(2:end);
end        

end
