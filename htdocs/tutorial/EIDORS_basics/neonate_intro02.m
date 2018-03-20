vv= eidors_readdata('P04P-1016.get'); vi=vv(:,45); vh=vv(:,61);
imr = inv_solve(imdl,vh,vi);

clf; axes('position',[0.05,0.5,0.25,0.45]);
imr.calc_colours = struct('ref_level',0,'greylev',0.2,'backgnd',[1,1,1]);
show_slices(imr);
