function [mdl] = reordernodes(mdl)
%REORDERNODES
%MDL = reordernodes(MDL)
%
%Function to put the internal nodes first in the node list, and
%subsequently renumber the nodes in mdl.elem(i) for each element

%Cache the node, element and boundary structure
boundstruc=mdl.boundary; nodestruc=mdl.nodes; 

%No. of nodes, elements, nodes per element and problem dimension
nnodes=size(nodestruc,1); nodedim=size(nodestruc,2);

%Use Unique to give all node no's on boundary and no. bound/int nodes
boundunique=unique(boundstruc(:));
nbound=size(boundunique,1); nint=nnodes-nbound;

%Initialise boudnary/internal arrays (with extra column for global numbers)
boundary=zeros(nbound,nodedim+1); internal=zeros(nint,nodedim+1);
boundtest=1; inttest=1; %Counters of boundary/internal numbers

for j=1:nnodes
    if(boundtest<=nbound)
        if(j==boundunique(boundtest))
            boundary(boundtest,:)=[j,nodestruc(j,:)]; %List boundary vtx numbers/coords
            boundtest=boundtest+1;
        else
            internal(inttest,:)=[j,nodestruc(j,:)]; %List internal vtc numbers/coords
            inttest=inttest+1;
        end
    else
        internal(inttest,:)=[j,nodestruc(j,:)]; %List internal vtc numbers/coords
        inttest=inttest+1;
    end
end 

%Store the list of boundary and internal nodes (with their global numbers)
mdl.solver.boundary=boundary; mdl.solver.internal=internal;
mdl.solver.int=nint; mdl.solver.bound=nbound;

% Change the ordering so that internal nodes first then boundary. Keep an
% index vector of the original node numbers
nodeidx=[mdl.solver.internal(:,1);mdl.solver.boundary(:,1)];
nodemat=[mdl.solver.internal(:,2:nodedim+1);mdl.solver.boundary(:,2:nodedim+1)];

%Change global node numbers in mdl.elem
%Cache the elem structure with a copy and find no. of elements/nodes
elemstruc=mdl.elem; elemstruccopy=elemstruc; 
nelems=size(elemstruc,2); nnodes=size(nodestruc,1);

for nn=1:nnodes %Loop over nodes
    ii=nodeidx(nn); %The original node number
    for i=1:nelems %Loop over elements
        nnodeselems=size(elemstruc(i).nodes,2); %No. nodes/element
        for j=1:nnodeselems
            if(elemstruccopy(i).nodes(j)==ii) %Test if old node no.
                elemstruc(i).nodes(j)=nn; %If so then swap for new number
            end
        end
    end
end

%Change global node numbers in mdl.boundary
%Cache boundary strucure and find no. of boundary edges(faces)/nodes
boundstruc=mdl.boundary; boundstruccopy=boundstruc;
nboundedge=size(boundstruc,1); nnodes=size(nodestruc,1);

for jj=1:nnodes %Loop over nodes
    qq=nodeidx(jj); %The original node number
    for b=1:nboundedge %Loop over boundary edges
        nnodesbound=size(boundstruc,2); %No. of nodes/boundary edge
        for c=1:nnodesbound
            if(boundstruccopy(b,c)==qq) %Test if old node no.
                boundstruc(b,c)=jj; %If so then swap for new number
            end
        end
    end
end

%Put changed element and node matrix back into mdl
mdl.elem=elemstruc; mdl.nodes=nodemat; mdl.boundary=boundstruc;

end