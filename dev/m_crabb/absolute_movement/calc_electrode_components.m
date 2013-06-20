function elec_comp = calc_electrode_components(fwd_model)

%This function calculates the following on each electrode, and stores these
%in the structure elec_comp{i}.NAME
%1. COM - The centre of mass of electrode
%2. NORMAL - The (average by area) normal vector to the electrode
%3. TANGENT - The (average by area) tangent space of the electrode
          %2D - Single vector - Defined so these are clockwise basis
          %3D - Two vectors - Defined for some basis. FIGURE THIS OUT!!!!
%4. BOUNDARY_NODES - List of electrode's boundary nodes
          %2D - Simply two nodes numbers on each side.
          %3D - All the nodes on the boundary DO I NEED CONNECTIONS HERE??
%5. BOUNDARY_NODES_DIR - Normal to electrode at each boundary node. We
%have tangent space basis of electrode, and so to reasonable level of 
%approximation, we can say this is linear combination of tangent space
%basis vectors i.e. matrix - n_boundary_nodes*(ndim-1)
          %2D - Length two vector corresponding to boundary_nodes with
          %values +/-1 corresponding to +/- tangent basis
          %3D - This is a matrix. Each column is for given boundary_nodes,
          %and we have \alpha, \beta i.e. vdE = \alpha t_{1} + \beta t_{2}
%
%TODO - 1. Calculate the normals in linear_bound_reorder
%     - 2. Gram-Schmidt orthogonalisation for tangent space?
%          
          
%Testing to chec orientation          
testing=1;

%Get the boundaries so that they are in a consistent numbering
fwd_model=linear_bound_reorder(fwd_model,-1);

%Get all the nodes, elems, elecs and boundary info
elecstruc =  fwd_model.electrode; n_elec = length(elecstruc);
boundstruc = fwd_model.boundary; n_bound = length(boundstruc(:,2));
nodestruc =  fwd_model.nodes; n_node = length(nodestruc(:,2)); nodedim = length(nodestruc(1,:));
elemstruc =  fwd_model.elems; n_elem = length(elemstruc(:,2));

%For each electrode
for i=1:n_elec
    %Calculate the centre of mass of the electrode   
    elec_nodes = elecstruc(i).nodes;
    pos = mean(nodestruc(elec_nodes,:),1);
    %Store the centre of mass in structure
    elec_comp{i}.com = pos;
    
    %Find which boundarys electrode has
    [bdy_idx,bdy_area]=find_electrode_bdy(boundstruc(:,1:nodedim), nodestruc,elecstruc(i).nodes);        
    boundidx_i=bdy_idx; area_i=bdy_area;
    n_bound_i = length(boundidx_i);     
    
    %Calculate the (area averaged) normal to boundary and tangent
    norm_i=zeros(nodedim,1);

    %Loop on electrode boundaries - area weighted electrode normal/tangent
    for j=1:n_bound_i
        %The coordinates of this boundarys nodes
        thisb=nodestruc(boundstruc(boundidx_i(j),:),:);        

        %Find its outer unit normal and its area  
        if(nodedim==2)            
            %Calculate the area and normalised clockwise tangent vector
            area_j = norm( thisb(1,:)-thisb(2,:) );           
            tang_j = (thisb(2,:)-thisb(1,:))'/area_j;   
            
            %Outward normal vector is 90 anticlockwise rotation
            norm_j = [0,-1;1,0]*tang_j;
        elseif(nodedim==3)           
            %Calculate the normalised outer normal to surface
            norm_j = cross(thisb(3,:)-thisb(1,:),thisb(2,:)-thisb(3,:))';
            area_j = norm(norm_j);                   
        end
        %Add on to normal weighted by its area
        norm_i = norm_i + norm_j*area_j;   
    end  
    
    %Normalise the normal
    norm_i = norm_i/norm(norm_i);  
    
    %Calculate the tangents to normal
    if(nodedim==2)
        %Average tangent vector is 90 rotation clockwise of normal
        tang_i = [0,1;-1,0]*norm_i;        
    elseif(nodedim==3)       
        %Form orthogonal basis on plane defined by normal - case by case
        %Gram-Schmidt procedure??? using 0.2 for some stability
        if(abs(norm_i(1))>0.2)
           tang_i1 = [-norm_i(2)/norm_i(1),1,0]; tang_i1=tang_i1/norm(tang_i1);
           tang_i2=cross(tang_i1,norm_i); tang_i2=tang_i2/norm(tang_i2);                      
        elseif(abs(norm_i(2))>0.2)
           tang_i1 = [-norm_i(2)/norm_i(1),1,0]; tang_i1=tang_i1/norm(tang_i1);
           tang_i2=cross(tang_i1,norm_i); tang_i2=tang_i2/norm(tang_i2);               
        elseif(abs(norm_i(3))>0.2)
           tang_i1 = [-norm_i(2)/norm_i(1),1,0]; tang_i1=tang_i1/norm(tang_i1);
           tang_i2=cross(tang_i1,norm_i); tang_i2=tang_i2/norm(tang_i2);                           
        else
            error('Can not compute normals correctly')
        end                
        
        %Now put this tangents as two columns of a matrix
        tang_i=[tang_i1;tang_i2]';
    end
    
    %Store the tangent(s) and normal in the structure
    elec_comp{i}.tangent = tang_i; %Tangent vectors by column
    elec_comp{i}.normal = norm_i; %Normal vector by column
    
    %Find the nodes on end points of electrode 
    if(nodedim==2)  
        %List all nodes on ith electrode as vector (not unique yet!)
        nodes_elec_i=[];
        for ii=1:length(boundidx_i)
            nodes_elec_i = [nodes_elec_i,boundstruc(boundidx_i(ii),:)];
        end
        %Find unique nodes and get a count histogram
        numbers=unique(nodes_elec_i); count=histc(nodes_elec_i,numbers); 
        
        %Loop through the unique nodes
        cnt=0;        
        for jj=1:length(count)        
           if(count(jj)==1)
              %If one occurence (i.e. it is on dEi) add to list
              cnt=cnt+1;               
              elec_comp{i}.boundary_nodes(cnt)=numbers(jj);
              %Get coordinates of the node
              nodecoord=nodestruc(numbers(jj),:);
              
              %Now compute t.(v-COM) to find vdE relative to t                                         
              alpha= elec_comp{i}.tangent'*(nodecoord-elec_comp{i}.com)';
              
              if(alpha>0)
                  %The multiple of the tangent vector is positive
                  elec_comp{i}.boundary_nodes_dir(cnt)=1; %alpha*t
              elseif(alpha<0)
                  %This multiple of the tangent vector is negative
                  elec_comp{i}.boundary_nodes_dir(cnt)=-1; %alpha*t
              end            
           end
        end       
    elseif(nodedim==3)
        %3D we need to calculate the \alpha and \beta i.e. dot product, to
        %find the components (in the given basis of tangent space) of the
        %movement
               
        cnt=0;
        %Loop through each electrode boundary and find if node in
        for iii=1:length(bdy_idx)
           %Get the nodes of this boundary
           biiinodes=boundstruc(bdy_idx(iii),:);
           
           %Find boundary which shares a node with this boundary
           [rn1,blah]=find(boundstruc==biiinodes(1));
           [rn2,blah]=find(boundstruc==biiinodes(2));
           [rn3,blah]=find(boundstruc==biiinodes(3));

           %Lists of boundaries that have two nodes in common (two!!)
           rn12=intersect(rn1,rn2); %node 1/2                      
           rn13=intersect(rn1,rn3); %node 1/3          
           rn23=intersect(rn2,rn3); %node 2/3                                                          
               
           %Set difference with electrode boundaries to find
           rn12diff = setdiff(rn12,bdy_idx);
           rn13diff = setdiff(rn13,bdy_idx);
           rn23diff = setdiff(rn23,bdy_idx);           
           
           %Go through cases if one of above not empty, store this
           %information inf elec_comp structure
           if(~isempty(rn12diff))
               cnt=cnt+1;
               elec_comp{i}.partialEi_bdyidx(cnt)=bdy_idx(iii);
               elec_comp{i}.partialEi_edge(cnt,1)=biiinodes(1);
               elec_comp{i}.partialEi_edge(cnt,2)=biiinodes(2);              
           elseif(~isempty(rn13diff))
               cnt=cnt+1;
               elec_comp{i}.partialEi_bdyidx(cnt)=bdy_idx(iii);
               elec_comp{i}.partialEi_edge(cnt,1)=biiinodes(1);
               elec_comp{i}.partialEi_edge(cnt,2)=biiinodes(3);                                 
           elseif(~isempty(rn23diff))
               cnt=cnt+1;
               elec_comp{i}.partialEi_bdyidx(cnt)=bdy_idx(iii);
               elec_comp{i}.partialEi_edge(cnt,1)=biiinodes(2);
               elec_comp{i}.partialEi_edge(cnt,2)=biiinodes(3);                  
           else
               %Do nothing - boundary is internal to electrode
           end                                       
        end                
    end        
end

end

function [mdl]=linear_bound_reorder(mdl,ccw)    
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