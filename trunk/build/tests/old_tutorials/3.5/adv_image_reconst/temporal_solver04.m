% Image reconstruction of moving objects $Id$

time_steps=  3; ts_expand= 5;
time_weight= .8;
ts_vec= -time_steps:time_steps;

image_select= .25*length(xyr_pt)+1:ts_expand:.75*length(xyr_pt)+1;

% GN Solver
 vi_sel = vi(:,image_select);
 img= inv_solve( imdl_GN, vh, vi_sel);
 animate_reconstructions('temporal_solver04a', img);

% Temporal Solver
 k=1;
 for i= image_select
   vi_sel= vi(:,i+ts_vec);
   imgs= inv_solve( imdl_TS, vh, vi_sel);
   img.elem_data(:,k)= imgs.elem_data(:,1+time_steps);
   k=k+1;
 end
 animate_reconstructions('temporal_solver04c', img);


% Kalman Solver
 vi_sel = vi(:,image_select);
 img= inv_solve( imdl_KS, vh, vi_sel);
 animate_reconstructions('temporal_solver04d', img);
