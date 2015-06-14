function [boundary,elems,nodes]=fem_1st_to_higher_order(fwd_model)
% FEM_1ST_TO_HIGH_ORDER:  Modify the FEM for high order FEM called as
%    [bound,elem,nodes] = fem_1st_to_higher_order( fwd_model )
% where field fwd_model.approx_type = 'tri3'  - 2D linear
%                                   = 'tri6'  - 2D quadratic
%                                   = 'tri10' - 2D cubic
%                                   = 'tet4'  - 3D linear
%                                   = 'tet10' - 3D quadratic
% fwd_model : is a fwd_model structure
% boundary/elems : new boundary and element numbers in connectivity matrix
% nodes : fem nodes including extra nodes added in refinement
%
%M Crabb - 29.06.2012

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return ; end

copt.cache_obj = {fwd_model.nodes, fwd_model.elems, ...
             fwd_model.approx_type, fwd_model.boundary};
copt.fstr = 'boundaryN_elemsN_nodes';

[boundary,elems,nodes]= eidors_cache(@mc_fem_modify,fwd_model,copt);
end

function [boundary,elems,nodes]=mc_fem_modify(mdl)
%Function  that takes a fwd_model structure (with at least fwd_model.boundary,
%fwd_model.elems, fwd_model.nodes), reorders these to form anti-clockwise sets, and
%then adds extra nodes for p-refinement.

%Reorder element nodes anti-clockwise relative to first node
mdl = linear_elem_reorder(mdl,-1);

%Reorder boundary nodes anti-clockwise relative to interior node
mdl = linear_bound_reorder(mdl,-1);

%Find nodes, elems, boundary and problem dimensions
n_nodes = size(mdl.nodes,1); n_elems = size(mdl.elems,1); n_bounds=size(mdl.boundary,1);
nodedim = size(mdl.nodes,2);

%Copy nodes, elems and boundaries
nodes_n=mdl.nodes; boundary_n=mdl.boundary; elems_n=mdl.elems;

%Compare strings
if(strcmp(mdl.approx_type,'tri6'))    
  %No nodes (really vertices) per element and boundary
  n_nodes_per_elem = 3; 
  n_nodes_per_boundary=1;
elseif(strcmp(mdl.approx_type,'tri10'))
  %No nodes (really vertices) per element and boundary
  n_nodes_per_elem = 7; 
  n_nodes_per_boundary=2;
elseif(strcmp(mdl.approx_type,'tet10'));
  %No nodes (really vertices) per element and boundary
  n_nodes_per_elem = 6; 
  n_nodes_per_boundary=3;
else
    error('mc_fem_modify: Element type ("%s") not recognised',mdl.approx_type);
end

%Maximum number of new nodes i.e. unconnected elemets and boundaries
n_nodes_e = n_nodes + n_nodes_per_elem*n_elems;
n_nodes_e_b = n_nodes_e + n_nodes_per_boundary*n_bounds;

%Now reshape the nodes, elems and boundary for extra FEM nodes
nodes_n(n_nodes+1:n_nodes_e_b,:) = 0;
elems_n(:,nodedim+2:nodedim+1+n_nodes_per_elem) = reshape((n_nodes+1):n_nodes_e,n_nodes_per_elem,n_elems)';
boundary_n(:,nodedim+1:nodedim+n_nodes_per_boundary) = reshape((n_nodes_e+1):n_nodes_e_b,n_nodes_per_boundary,n_bounds)';

%Calculate positions of now nodes dependent on approx_type
if(strcmp(mdl.approx_type,'tri6'))    
  nodes_n(elems_n(:,4),:)  = 0.5*(nodes_n(elems_n(:,1),:) + nodes_n(elems_n(:,2),:));
  nodes_n(elems_n(:,5),:)  = 0.5*(nodes_n(elems_n(:,2),:) + nodes_n(elems_n(:,3),:));
  nodes_n(elems_n(:,6),:)  = 0.5*(nodes_n(elems_n(:,3),:) + nodes_n(elems_n(:,1),:));
  nodes_n(boundary_n(:,3),:) = 0.5*(nodes_n(boundary_n(:,1),:) + nodes_n(boundary_n(:,2),:));
elseif(strcmp(mdl.approx_type,'tri10'))
  nodes_n(elems_n(:,4),:) = 2.0/3.0*nodes_n(elems_n(:,1),:) + 1.0/3.0*nodes_n(elems_n(:,2),:);
  nodes_n(elems_n(:,5),:) = 1.0/3.0*nodes_n(elems_n(:,1),:) + 2.0/3.0*nodes_n(elems_n(:,2),:);
  nodes_n(elems_n(:,6),:) = 2.0/3.0*nodes_n(elems_n(:,2),:) + 1.0/3.0*nodes_n(elems_n(:,3),:);
  nodes_n(elems_n(:,7),:) = 1.0/3.0*nodes_n(elems_n(:,2),:) + 2.0/3.0*nodes_n(elems_n(:,3),:);
  nodes_n(elems_n(:,8),:) = 2.0/3.0*nodes_n(elems_n(:,3),:) + 1.0/3.0*nodes_n(elems_n(:,1),:);
  nodes_n(elems_n(:,9),:) = 1.0/3.0*nodes_n(elems_n(:,3),:) + 2.0/3.0*nodes_n(elems_n(:,1),:);
  nodes_n(elems_n(:,10),:) = 1.0/3.0*nodes_n(elems_n(:,1),:) + 1.0/3.0*nodes_n(elems_n(:,2),:) + 1.0/3.0*nodes_n(elems_n(:,3),:);
  nodes_n(boundary_n(:,3),:) = 2.0/3.0*nodes_n(boundary_n(:,1),:) + 1.0/3.0*nodes_n(boundary_n(:,2),:);
  nodes_n(boundary_n(:,4),:) = 1.0/3.0*nodes_n(boundary_n(:,1),:) + 2.0/3.0*nodes_n(boundary_n(:,2),:);
elseif(strcmp(mdl.approx_type,'tet10'))
  nodes_n(elems_n(:,5),:) = 0.5*nodes_n(elems_n(:,1),:) + 0.5*nodes_n(elems_n(:,2),:);
  nodes_n(elems_n(:,6),:) = 0.5*nodes_n(elems_n(:,1),:) + 0.5*nodes_n(elems_n(:,3),:);
  nodes_n(elems_n(:,7),:) = 0.5*nodes_n(elems_n(:,1),:) + 0.5*nodes_n(elems_n(:,4),:);
  nodes_n(elems_n(:,8),:) = 0.5*nodes_n(elems_n(:,2),:) + 0.5*nodes_n(elems_n(:,3),:);
  nodes_n(elems_n(:,9),:) = 0.5*nodes_n(elems_n(:,3),:) + 0.5*nodes_n(elems_n(:,4),:);
  nodes_n(elems_n(:,10),:) = 0.5*nodes_n(elems_n(:,2),:) + 0.5*nodes_n(elems_n(:,4),:);
  nodes_n(boundary_n(:,4),:) = 0.5*nodes_n(boundary_n(:,1),:) + 0.5*nodes_n(boundary_n(:,2),:);
  nodes_n(boundary_n(:,5),:) = 0.5*nodes_n(boundary_n(:,1),:) + 0.5*nodes_n(boundary_n(:,3),:);
  nodes_n(boundary_n(:,6),:) = 0.5*nodes_n(boundary_n(:,2),:) + 0.5*nodes_n(boundary_n(:,3),:);
else
    error('mc_fem_modify: Element type ("%s") not recognised',mdl.approx_type);
end

%Find the unique nodes by row 
[nodes, a, b] = unique(nodes_n(n_nodes+1:end,:),'rows');

%Add the unique nodes to the model
nodes_n = [mdl.nodes; nodes];

%Now find unique nodes
c = [1:n_nodes n_nodes+b']; 
elems_n = c(elems_n);
boundary_n = c(boundary_n);

%Reassign matrices
nodes=nodes_n; elems=elems_n; boundary=boundary_n;

end

%Function to reorder nodes consistently within element/boundary
%Elem Default : Arrange nodes locally in mdl.elems so that:
%2D: Counter-clockwise looking down onto element
%3D: Counter-clockwise looking into element from first node
%    (clockwise looking down from ouside element opposite first node)
%Boundary Default: Arrange nodes locally in mdl.boundary so that:
%2D: Counter-clockwise looking down relative to interior node
%3D: Clockwise relative looking into element from interior node
%    (counter-clockwise looking down from outside element opposite interior node)
%Conventions MUST match reference element geometry in system_mat_higher_order
function [mdl] = linear_elem_reorder(mdl,ccw)
    %Find elements, nodes and no. of elements
    elemstruc=mdl.elems;  nodestruc=mdl.nodes; nelems=size(elemstruc,1);
    for e=1:nelems;
        %Row vector of the node numbers and the no. of vertices
        enodes = elemstruc(e,:); elenode = size(enodes,2);
    
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
        elemstruc(e,:)=enodes; %Put enodes back into elementnodes matrix
    end
    %Reassign the elements
    mdl.elems=elemstruc;
end

function [mdl]=linear_bound_reorder(mdl,ccw)
    %Reorder nodes on each electrode's boundary
    %SPEED UP
    %1fix_model with option opt.bound_elec2elem=1
    
    %Find boundary, elems, nodes, elecs
    boundstruc=mdl.boundary; elemstruc=mdl.elems; nodestruc=mdl.nodes; elecstruc=mdl.electrode;
    %Find no. elecs and node dimension
    nelecs=size(elecstruc,2); nodedim=size(nodestruc,2);

    %Reorder boundaries belonging to electrodes
    for ke=1:nelecs
        %The boundary numbers and areas, outputs rows of mdl.boundary of electrode
        bdy_idx=find_electrode_bdy(boundstruc(:,1:nodedim),nodestruc,elecstruc(ke).nodes);                 
            
        for ii=1:length(bdy_idx);  
            %Get bdy idx
            bb=bdy_idx(ii);
                                    
            %Vector of vertes numbers of this boundary
            bbnodes=boundstruc(bb,:);
            if(nodedim==2) %2D problem
                %Find row(s) of elems for which each boundary vertex belongs
                [rownode1,blah]=find(elemstruc==bbnodes(1));
                [rownode2,blah]=find(elemstruc==bbnodes(2));
        
                %Intersection of rownode1 and rownode2 is the unique element
                boundiielem=intersect(rownode1,rownode2);
            elseif(nodedim==3) %3D problem
                %Find row(s) of elems for which each node belongs
                [rownode1,blah]=find(elemstruc==bbnodes(1));
                [rownode2,blah]=find(elemstruc==bbnodes(2));
                [rownode3,blah]=find(elemstruc==bbnodes(3));

                %Intersection of rownode1 and rownode2 gives a choic
                rownode1node2=intersect(rownode1,rownode2);

                %Intersection of rownode3 with vector above is unique element
                boundiielem=intersect(rownode3,rownode1node2);  
            end
            %Store this unique number in a vecto
            %FIXME!!! 3D Inclusion models seem to include the boundaries of the
            %inclusions as boundaries. In this case there may be multiple
            %boundaries. Just choose first one for time being, and if this
            %is ACTUAL boundary, then should be unique.....
            beleind=boundiielem(1);     
            
            %Coordinate of nodes of boundaries element
            belenodes=elemstruc(beleind,:);
  
            %Find unique internal node and order element node so internal is first 
            intnode=setdiff(belenodes,bbnodes); elemboundnode=[intnode,bbnodes];
    
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
            boundstruc(bb,:)=elemboundnode(2:end);
        end              
    end
    %Reassign the boundary
    mdl.boundary=boundstruc;
end

function do_unit_test
     do_unit_test_2D;
     do_unit_test_3D;
end

function do_unit_test_2D
    imdl=mk_common_model('c2C',16);
    fmdl=imdl.fwd_model;
    fmdl.approx_type='tri6';
    [bou,ele,nod]=fem_1st_to_higher_order(fmdl);

    unit_test_cmp('NOD',nod(1:5,:), ...
[                  0                   0
  -0.005450260769179   0.083154910269884
   0.083154910269884   0.005450260769179
   0.005450260769179  -0.083154910269884
  -0.083154910269884  -0.005450260769179],1e-8);
    unit_test_cmp('ELE',ele(1:5,1:6), ...
[1, 3, 2, 783, 780, 755
 1, 4, 3, 760, 785, 783
 1, 5, 4, 732, 735, 760
 1, 2, 5, 755, 730, 732
 7, 2, 3, 792, 780, 810]);

end

function do_unit_test_3D
    fmdl=mk_library_model('adult_male_16el');
    fmdl.stimulation = mk_stim_patterns(16,1,'{ad}','{ad}');
    fmdl.approx_type='tet10';
    [bou,ele,nod]=fem_1st_to_higher_order(fmdl);

    unit_test_cmp('NOD',nod(1:5,:), ...
[ -0.812999999999999, -0.439800000000001, 0
  -0.720700000000001, -0.494600000000000, 0
  -0.886900000000000, -0.361400000000000, 0
  -0.616600000000001, -0.523199999999999, 0
  -0.516400000000002, -0.562899999999999, 0],1e-8 );

end
