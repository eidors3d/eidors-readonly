function [img]=mc_calc_system_mat_CEM(img)
%Assemble the total stiffness matrix

%Get the forward model
mdl=img.fwd_model;

%Find electrode stucture and no.of electrodes and initialize vector
elecstruc=mdl.electrode; nelecs=size(elecstruc,2);

%Find node structure and find no.nodes 
nodestruc=mdl.nodes; nnodes=size(nodestruc,1); 

%Assemble the total system matrix
Am=mdl.solver.Am; Az=mdl.solver.Az; Aw=mdl.solver.Aw; Ad=mdl.solver.Ad;

%Total system matrix: At = ( Am+Az , Aw )
%                          (  Aw'  , Ad )
At=zeros(nnodes+nelecs,nnodes+nelecs);

At(1:nnodes,1:nnodes) = Am+Az;
At(1:nnodes,nnodes+1:nnodes+nelecs) = Aw;
At(nnodes+1:nnodes+nelecs,1:nnodes)=Aw';
At(nnodes+1:nnodes+nelecs,nnodes+1:nnodes+nelecs)=Ad;

%Put the total system matrix in forward model and put this into image
mdl.solver.At=At; img.fwd_model=mdl;