% Data: eidors3d.sf.net/data_contrib/if-neonate-spontaneous/if-neonate-spontaneous.zip
vv= eidors_readdata('P04P-1016.get');

% solve with reference to the mean
imgall = inv_solve(imdl,mean(vv,2),vv);
[insp, expi] = find_frc(imgall,[],13,[],2); % find breaths

% use expirations as reference
imgr = inv_solve(imdl,vv(:,expi),vv(:,insp(2:end)));
   imgr.calc_colours.ref_level= 0;
   imgr.calc_colours.backgnd= [1 1 1];
   imgr.calc_colours.greylev= 0.001;
   imgr.show_slices.img_cols = 4;
   imgr.show_slices.sep      = 2;
show_slices(imgr);

print_convert neonate_ex03a.png
