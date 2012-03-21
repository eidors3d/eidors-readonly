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