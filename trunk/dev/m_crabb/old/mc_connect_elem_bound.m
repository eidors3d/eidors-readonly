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
    
        if(nodedim==2) %2D problem
            %Find row(s) of elems for which each boundary vertex belongs
            [rownode1,~]=find(elemstrucold==boundaryiinodes(1));
            [rownode2,~]=find(elemstrucold==boundaryiinodes(2));
        
            %Intersection of rownode1 and rownode2 is the unique element
            boundiielem=intersect(rownode1,rownode2);
        elseif(nodedim==3) %3D problem
            %Find row(s) of elems for which each node belongs
            [rownode1,~]=find(elemstrucold==boundaryiinodes(1));
            [rownode2,~]=find(elemstrucold==boundaryiinodes(2));
            [rownode3,~]=find(elemstrucold==boundaryiinodes(3));

            %Intersection of rownode1 and rownode2 gives a choice
            rownode1node2=intersect(rownode1,rownode2);

            %Intersection of rownode3 with vector above is unique element
            boundiielem=intersect(rownode3,rownode1node2);         
        end
        %Store this unique number in a vector
        bound_elem(ii)=boundiielem;
    end    
end