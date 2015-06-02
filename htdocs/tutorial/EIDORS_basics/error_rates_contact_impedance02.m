%Loop over contact impedances
z_c=[0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000];
for z_cj=1:length(z_c) % z_c=[0.1,1000];
%Reinitialise model and make image
mdl=mdl_init; img=mk_image(mdl,1);
mdl.electrode(1).z_contact = z_c(z_cj); 
img.fwd_model.electrode(1).z_contact=z_c(z_cj);
mdl.electrode(2).z_contact = z_c(z_cj); 
img.fwd_model.electrode(2).z_contact=z_c(z_cj);

%Solve forward problem with linear, quadratic, cubic and find H1-error
tic;
  fprintf(1,'H1 Lin');
  [lin_errL2(1,z_cj),lin_errH1(1,z_cj),lin_errH1tot(1,z_cj),lin_errI(1,z_cj), ...
   lin_errUS(1,z_cj),lin_errUM(1,z_cj),lin_errUSM(1,z_cj), ...
   lin_timing_solver(1,z_cj),lin_DOF(1,z_cj)]= error_2D_squ_CEM(img,'tri3',1);
toc;
tic;
  fprintf(1,'H1 Quad');
  [quad_errL2(1,z_cj),quad_errH1(1,z_cj),quad_errH1tot(1,z_cj),quad_errI(1,z_cj), ...
  quad_errUS(1,z_cj),quad_errUM(1,z_cj),quad_errUSM(1,z_cj), ...
  quad_timing_solver(1,z_cj),quad_DOF(1,z_cj)]= error_2D_squ_CEM(img,'tri6',1);
toc;
tic; fprintf(1,'H1 Cub'); [cub_errL2(1,z_cj),cub_errH1(1,z_cj),cub_errH1tot(1,z_cj),cub_errI(1,z_cj),cub_errUS(1,z_cj),cub_errUM(1,z_cj),cub_errUSM(1,z_cj),cub_timing_solver(1,z_cj),cub_DOF(1,z_cj)]=error_2D_squ_CEM(img,'tri10',1); toc;           

%Loop over different refinements
for ii=1:nmesh-1
    %Refine model uniformly and find errors with lin, quad, cubic                
    mdl=h_refine(mdl);  img=mk_image(mdl,1);        
    tic; fprintf(1,'H1 Lin'); [lin_errL2(ii+1,z_cj),lin_errH1(ii+1,z_cj),lin_errH1semi(ii+1,z_cj),lin_errI(ii+1,z_cj),lin_errUS(ii+1,z_cj),lin_errUM(ii+1,z_cj),lin_errUSM(ii+1,z_cj),lin_timing_solver(ii+1,z_cj),lin_DOF(ii+1,z_cj)]=error_2D_squ_CEM(img,'tri3',0); toc;
    tic; fprintf(1,'H1 Quad'); [quad_errL2(ii+1,z_cj),quad_errH1(ii+1,z_cj),quad_errH1semi(ii+1,z_cj),quad_errI(ii+1,z_cj),quad_errUS(ii+1,z_cj),quad_errUM(ii+1,z_cj),quad_errUSM(ii+1,z_cj),quad_timing_solver(ii+1,z_cj),quad_DOF(ii+1,z_cj)]=error_2D_squ_CEM(img,'tri6',0); toc;
    tic; fprintf(1,'H1 Cub'); [cub_errL2(ii+1,z_cj),cub_errH1(ii+1,z_cj),cub_errH1semi(ii+1,z_cj),cub_errI(ii+1,z_cj),cub_errUS(ii+1,z_cj),cub_errUM(ii+1,z_cj),cub_errUSM(ii+1,z_cj),cub_timing_solver(ii+1,z_cj),cub_DOF(ii+1,z_cj)]=error_2D_squ_CEM(img,'tri10',0); toc;           
end

%Refine model uniformly and find errors with lin               
mdl=h_refine(mdl);  img=mk_image(mdl,1);        
tic; fprintf(1,'H1 Lin'); [lin_errL2(nmesh+1,z_cj),lin_errH1(nmesh+1,z_cj),lin_errH1semi(nmesh+1,z_cj),lin_errI(nmesh+1,z_cj),lin_errUS(nmesh+1,z_cj),lin_errUM(nmesh+1,z_cj),lin_errUSM(nmesh+1,z_cj),lin_timing_solver(nmesh+1,z_cj),lin_DOF(ii+1,z_cj)]=error_2D_squ_CEM(img,'tri3',0); toc;       
end
