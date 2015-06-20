vv= eidors_readdata('p07b.get');

% solve with reference to the mean
imgall = inv_solve(imdl,mean(vv,2),vv);
[insp, expi] = find_frc(imgall,[],13,[],2); % find breaths

% use expirations as reference
imgr = inv_solve(imdl,vv(:,expi),vv(:,insp));
   imgr.calc_colours.ref_level= 0;
   imgr.calc_colours.backgnd= [1 1 1];
   imgr.calc_colours.greylev= 0.001;
   imgr.show_slices.img_cols = 4;
   imgr.show_slices.sep      = 2;
   imgr.calc_colours.cmap_type = 'swisstom';
show_slices(imgr);

print_convert neonate_ex03a.png
