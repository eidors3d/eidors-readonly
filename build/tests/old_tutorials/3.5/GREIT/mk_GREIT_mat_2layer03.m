% $Id$

rimg = inv_solve(imdl,vh,vi); %Reconstruct
rimg.calc_colours.ref_level=0;

show_fem(rimg); axis square;

print_convert mk_GREIT_mat_2layer03a.png '-density 60'
