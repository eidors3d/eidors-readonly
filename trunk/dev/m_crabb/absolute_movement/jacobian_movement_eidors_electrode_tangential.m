function [J] = jacobian_movement_eidors_electrode_tangential(fwd_model,img)
%Now we want to loop through the electrode and find out the unique normal
%components
%1. For each electrode find its surface triangle
%2. For each triangle work out its outer unit normal, the normal of the
%electrode should then by the \sum_{i} n_{i} A_{i} / \sum_{i} A_{i}
%3. We then have the unique weighted normal
%4. find the tangenetial plane to this choose a direction along this 
%    2D - this is a line and we just choose a direction (maybe we can do
%    this so it is a right handed set?
%    3D - this is a plane, maybe make a consistent orientaiton
%
%In inver solver we somehow want to store all this tangential and normal
%infgormation on a field, so we do not need to work it out in a script

%Get the Jacobian matrix
J=jacobian_movement(img);

%TODO We want to reorder the boundary nodes so that they form a consistent set
% img.fwd_model=linear_bound_reorder(img.fwd_model);

%Get all the nodes, elems, elecs and boundary info
elecstruc = img.fwd_model.electrode; n_elec = length(elecstruc);
elemstruc = img.fwd_model.elems;     n_elem=length(elemstruc);
nodestruc =  img.fwd_model.nodes;    nodedim = length(nodestruc(1,:));

%Calculate the normal and tangential space of each electrode
elec_comp = calc_electrode_components(img.fwd_model);

%Ok we have split into tangential and normal components
J_move = J(:,n_elem+1:end);

for i=1:n_elec
    if(nodedim==2)
        %Now we get the two columns for this electrode i.e. in x and y
        xd = J_move(:,i);
        yd = J_move(:,i+n_elec);
        xyd=[xd,yd];
    
        %Now get the tangential part
        td = xyd*elec_comp{i}.tangent;
    
        %Put these into a movement for noraml and tangential componeonents
        J_move_comp(:,i)=td;            
    elseif(nodedim==3)
        %Now we get the two columns for this electrode i.e. in x and y
        xd = J_move(:,i);
        yd = J_move(:,i+n_elec);
        zd = J_move(:,i+2*n_elec);
        xyzd=[xd,yd,zd];
    
        %Now get the tangential part
        t1d = xyzd*elec_comp{i}.tangent(:,1);
        t2d = xyzd*elec_comp{i}.tangent(:,2);
    
        %Put these into a movement for noraml and tangential componeonents
        J_move_comp(:,i)=t1d;     
        J_move_comp(:,i+n_elec)=t2d;                    
    end
end

%Now insert the componented into the original Jacobian
J(:,n_elem+1:n_elem+(nodedim-1)*n_elec)=J_move_comp;
J(:,n_elem+(nodedim-1)*n_elec+1:end)=[];

end