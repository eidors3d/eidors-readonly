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