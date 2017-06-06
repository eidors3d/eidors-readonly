function mdl = h_refine(mdl)
%Perform uniform h refinement of square CEM

%Copy nodes, elems and boundaries
nodes=mdl.nodes; boundary=mdl.boundary; elems=mdl.elems;

%Find nodes, elems, boundary and problem dimensions
n_nodes = size(nodes,1); n_elems = size(elems,1); n_boundary=size(boundary,1);

%Maximum number of new nodes i.e. unconnected elemets and boundaries
n_nodes_new = n_nodes + 3*n_elems + n_boundary; %Maximum number of new nodes
n_boundary_new = n_boundary*2; %Number of new boundaries EXACT
n_elems_new = n_elems*4; %Number of new elements EXACT

%Storage for these new arrays
nodes_new=zeros(n_nodes_new,2);
elems_new=zeros(n_elems_new,3);
boundary_new=zeros(n_boundary_new,2);

%Keep original nodes in n_nodes_new
nodes_new(1:n_nodes,:)=nodes;

%Loop through each element and fine edge midpoints.
cnt=0;
for i=1:n_elems
%Find new nodes on midpoints
nodes_new(n_nodes+(i-1)*3+1,:) = 0.5*(nodes(elems(i,1),:) + nodes(elems(i,2),:));
nodes_new(n_nodes+(i-1)*3+2,:) = 0.5*(nodes(elems(i,1),:) + nodes(elems(i,3),:));
nodes_new(n_nodes+(i-1)*3+3,:) = 0.5*(nodes(elems(i,2),:) + nodes(elems(i,3),:));

%Create the four new elements
elems_new((i-1)*4+1,1) = elems(i,1); elems_new((i-1)*4+1,2) = n_nodes+3*(i-1)+1; elems_new((i-1)*4+1,3)=n_nodes+3*(i-1)+2;
elems_new((i-1)*4+2,1) = elems(i,2); elems_new((i-1)*4+2,2) = n_nodes+3*(i-1)+1; elems_new((i-1)*4+2,3)=n_nodes+3*(i-1)+3;
elems_new((i-1)*4+3,1) = elems(i,3); elems_new((i-1)*4+3,2) = n_nodes+3*(i-1)+2; elems_new((i-1)*4+3,3)=n_nodes+3*(i-1)+3;
elems_new((i-1)*4+4,1) = n_nodes+3*(i-1)+1; elems_new((i-1)*4+4,2) = n_nodes+3*(i-1)+2; elems_new((i-1)*4+4,3)= n_nodes+3*(i-1)+3;
end

%Loop through each boundary and fine edge midpoints
for i=1:n_boundary
%Find node and add a row of         
nodes_new(n_nodes+3*n_elems+i,:) = 0.5*(nodes(boundary(i,1),:) + nodes(boundary(i,2),:));

%Create the two new boundarys
boundary_new((i-1)*2+1,1) = boundary(i,1); boundary_new((i-1)*2+1,2) = n_nodes+3*n_elems+i; 
boundary_new((i-1)*2+2,1) = n_nodes+3*n_elems+i; boundary_new((i-1)*2+2,2) = boundary(i,2); 
end

%Find the unique nodes by row 
[nodes, a, b] = unique(nodes_new(n_nodes+1:end,:),'rows');

%Add the unique nodes to the model
nodes_new = [mdl.nodes; nodes];

%Now find unique nodes
c = [1:n_nodes n_nodes+b']; 
elems_new = c(elems_new);
boundary_new = c(boundary_new);

%Reassign matrices
nodes=nodes_new; elems=elems_new; boundary=boundary_new;

%Find new CEM electrode nodes
%Recalculate the number of nodes
n_nodes = size(nodes,1);
cntE1 = 0; cntE2 = 0;
%Loop through boundary and see if electrode is in between electrodes
for i=1:n_nodes
    if(nodes(i,2)==0) %We are on bottom boundary
        if( pi/5 <= nodes(i,1) && nodes(i,1)  <= 2*pi/5)
            cntE1 = cntE1+1;
            E1node(cntE1) = i;
        elseif(  3*pi/5 <= nodes(i,1) && nodes(i,1)  <= 4*pi/5)
            cntE2 = cntE2+1;
            E2node(cntE2) = i;                
        end
    end            
end
%Electrode
mdl.electrode(1).nodes=E1node;
mdl.electrode(2).nodes=E2node;        

%Reassign 
mdl.nodes=nodes; mdl.elems=elems; mdl.boundary=boundary;
    
end
