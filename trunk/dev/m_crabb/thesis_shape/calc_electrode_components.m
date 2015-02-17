function elec_comp = calc_electrode_components(fwd_model)
%Calculates the following storing in elec_comp{i}.NAME
%1. COM - Centre of mass of electrode
%2. NORMAL - The (average by area) normal vector to the electrode
%3. TANGENT - The (average by area) tangent space of the electrode
          %2D - Single vector - Defined so these are clockwise basis
          %3D - Two vectors - Defined for some basis.    
          
%Get the boundaries so that they are in a consistent numbering
fwd_model=linear_bound_reorder(fwd_model,-1);

%Get all the nodes, elems, elecs and boundary info
elecstruc =  fwd_model.electrode; n_elec = length(elecstruc);
boundstruc = fwd_model.boundary; 
nodestruc =  fwd_model.nodes; nodedim = length(nodestruc(1,:));

%For each electrode
for i=1:n_elec
    %Calculate the centre of mass of the electrode   
    elec_nodes = elecstruc(i).nodes;
    pos = mean(nodestruc(elec_nodes,:),1);
    %Store the centre of mass in structure
    elec_comp{i}.com = pos;
    
    %Find which boundarys electrode has
    [bdy_idx,bdy_area]=find_electrode_bdy...
        (boundstruc(:,1:nodedim), nodestruc,elecstruc(i).nodes);        
    boundidx_i=bdy_idx;
    n_bound_i = length(boundidx_i);     
    
    %Calculate the (area averaged) normal to boundary and tangent
    norm_i=zeros(nodedim,1);

    %Loop on elec bdys - area weighted normal/tangent
    for j=1:n_bound_i
        %The coordinates of this boundarys nodes
        thisb=nodestruc(boundstruc(boundidx_i(j),:),:);        

        %Find its outer unit normal and its area  
        if(nodedim==2)            
            %Calculate area and normalised clockwise tangent
            area_j = norm( thisb(1,:)-thisb(2,:) );           
            tang_j = (thisb(2,:)-thisb(1,:))'/area_j;   
            
            %Outward normal vector is 90 anticlockwise rotation
            norm_j = [0,-1;1,0]*tang_j;
        elseif(nodedim==3)           
            %Calculate the normalised outer normal to surface
            norm_j = cross(thisb(3,:)-thisb(1,:), ...
                thisb(2,:)-thisb(3,:))';
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
        %Case by case ortohognal basis on electrode
        if(abs(norm_i(1))>0.2)
           tang_i1 = [-norm_i(2)/norm_i(1),1,0]; 
           tang_i1=tang_i1/norm(tang_i1);
           tang_i2=cross(tang_i1,norm_i); 
           tang_i2=tang_i2/norm(tang_i2);                      
        elseif(abs(norm_i(2))>0.2)
           tang_i1 = [-norm_i(2)/norm_i(1),1,0]; 
           tang_i1=tang_i1/norm(tang_i1);
           tang_i2=cross(tang_i1,norm_i); 
           tang_i2=tang_i2/norm(tang_i2);               
        elseif(abs(norm_i(3))>0.2)
           tang_i1 = [-norm_i(2)/norm_i(1),1,0]; 
           tang_i1=tang_i1/norm(tang_i1);
           tang_i2=cross(tang_i1,norm_i); 
           tang_i2=tang_i2/norm(tang_i2);                           
        else
            error('Can not compute normals correctly')
        end                        
        %Now put this tangents as two columns of a matrix
        tang_i=[tang_i1;tang_i2]';
    end
    
    %Store the tangent(s) and normal in the structure
    elec_comp{i}.tangent = tang_i; %Tangent vectors by column
    elec_comp{i}.normal = norm_i; %Normal vector by column    
end
end

function [mdl]=linear_bound_reorder(mdl,ccw)    
    %Find boundary, elems, nodes, elecs
    boundstruc=mdl.boundary; elemstruc=mdl.elems; 
    nodestruc=mdl.nodes; elecstruc=mdl.electrode;
    %Find no. elecs and node dimension
    nelecs=size(elecstruc,2); nodedim=size(nodestruc,2);

    %Reorder boundaries belonging to electrodes
    for ke=1:nelecs
        %Boundary numbers/areas, output rows
        bdy_idx=find_electrode_bdy(boundstruc(:,1:nodedim), ...
            nodestruc,elecstruc(ke).nodes);                 
            
        for ii=1:length(bdy_idx);  
            %Get bdy idx
            bb=bdy_idx(ii);
                                    
            %Vector of vertes numbers of this boundary
            bbnodes=boundstruc(bb,:);
            if(nodedim==2) %2D problem
                %Row(s) of elems which each bdy vertex belongs
                [rownode1,blah]=find(elemstruc==bbnodes(1));
                [rownode2,blah]=find(elemstruc==bbnodes(2));
        
                %Intersectionr rownode1/rownode2 
                boundiielem=intersect(rownode1,rownode2);
            elseif(nodedim==3) %3D problem
                %Row(s) of elems which each bdy vertex belongs
                [rownode1,blah]=find(elemstruc==bbnodes(1));
                [rownode2,blah]=find(elemstruc==bbnodes(2));
                [rownode3,blah]=find(elemstruc==bbnodes(3));

                %Intersectionr rownode1/rownode2 
                rownode1node2=intersect(rownode1,rownode2);

                %Intersection rownode3 is unique element
                boundiielem=intersect(rownode3,rownode1node2);  
            end
            %Store this unique number in a vector
            beleind=boundiielem(1);     
            
            %Coordinate of nodes of boundaries element
            belenodes=elemstruc(beleind,:);
  
            %Unique internal nodes
            intnode=setdiff(belenodes,bbnodes); 
            elemboundnode=[intnode,bbnodes];
    
            %List by row coordinates of the element
            nodepos=nodestruc(elemboundnode,:); 
            nvertices=size(belenodes,2);     

            %Calculate area(2D)/volume(3D) from vertices
            area= det([ones(nvertices,1),nodepos]); 
            areasign=sign(area);
    
            %If positive area, swap the last two nodes
            if(areasign == -ccw) %Swap last two entries of enodes 
                temp=elemboundnode(end-1);
                elemboundnode(end-1)=elemboundnode(end);
                elemboundnode(end) = temp;
            end
               
            %Put the node numbers back into mdl.bound(bb) 
            boundstruc(bb,:)=elemboundnode(2:end);
        end              
    end
    %Reassign the boundary
    mdl.boundary=boundstruc;
end