% Test formation of complete electrode model

imdl= mk_common_model('a2c2',8);
fmdl= imdl.fwd_model;
for i=1:8; 
   fmdl.electrode(i).nodes = fmdl.electrode(i).nodes + [0,1];
end
img= calc_jacobian_bkgnd(imdl);

system_mat_1st_order(fmdl,img);
disp('success'); % if we have no errors
