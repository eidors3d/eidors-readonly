function [J_new] = jac_move_eidors_elec_tang_only(fwd_model,img)
%Calculate movement Jacobian projected onto normal and tangential
%components of each electrode

%Calculate movement Jacobian through eidors
J=jacobian_movement(img);

%Get all the nodes, elems, elecs and boundary info
elecstruc = img.fwd_model.electrode; n_elec = length(elecstruc);
elemstruc = img.fwd_model.elems;    
nodestruc =  img.fwd_model.nodes;    nodedim = length(nodestruc(1,:));

%Calculate the normal and tangential space of each electrode
elec_comp = calc_electrode_components(img.fwd_model);

%Ok we have split into tangential and normal components
J_move = J(:,end-3*n_elec+1:end); %c2f

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
        xd = J_move(:,i); yd = J_move(:,i+n_elec);  zd = J_move(:,i+2*n_elec);
        xyzd=[xd,yd,zd];
    
        %Now get the tangential part
        t1d = xyzd*elec_comp{i}.tangent(:,1);
        t2d = xyzd*elec_comp{i}.tangent(:,2);
    
        %Put these into a movement for noraml and tangential componeonents
        J_move_comp(:,i)=t1d;     
        J_move_comp(:,i+n_elec)=t2d;                    
    end
end

%We only want to keep the movement part
J_new=J_move_comp;
end