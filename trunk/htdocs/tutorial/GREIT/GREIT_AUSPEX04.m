s=1;
fnames = {'ReconstrMatrixGREITc.mat','ReconstrMatrixGREITt.mat'};
for i=1:length(fnames);
   load(fnames{i});
   imdl = mk_common_gridmdl('b2d',ReconstrMatrix');

   imgr= inv_solve(imdl, vh, vi);
   subplot(2,2,s+0); show_fem(imgr);
   title('Test Object Reconstruction');
   
   imgr.elem_data = mean(ReconstrMatrix==0)';
   subplot(2,2,s+1); show_fem(imgr);
   title('Domain Boundary');
s=s+2;end
print_convert GREIT_AUSPEX04a.png
