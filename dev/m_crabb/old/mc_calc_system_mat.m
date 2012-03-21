function [s_mat]=mc_calc_system_mat(fwd_model,img)
%Assemble the total stiffness matrix

%Find no. of electrodes and no. of ndoes
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 

%Test - Point/Complete electrodes. Assume no mixed model so test first elect
if(size(elecstruc(1).nodes,2)==1 && size(elecstruc(1).nodes,1)==1) %POINT ELECTRODE
    %IF POINT ELECTRODE
    [Am]=mc_calc_stiffness(fwd_model,img);
    s_mat=Am;
else %COMPLETE ELECTRODE
    [Am]=mc_calc_stiffness(fwd_model,img);
    [Aw,Az,Ad]=mc_calc_complete(fwd_model,img);
    s_mat=zeros(nnodes+nelecs,nnodes+nelecs);
    s_mat(1:nnodes,1:nnodes) = Am+Az;
    s_mat(1:nnodes,nnodes+1:nnodes+nelecs) = Aw;
    s_mat(nnodes+1:nnodes+nelecs,1:nnodes)=Aw';
    s_mat(nnodes+1:nnodes+nelecs,nnodes+1:nnodes+nelecs)=Ad;
end