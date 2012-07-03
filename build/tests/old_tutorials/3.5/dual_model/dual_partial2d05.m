% Dual Partial $Id$

clf; levels= [inf,inf,0,1,1];

% reconstruct fine model params onto coarse model
imgf= inv_solve(frec_mdl, vh, vi);
imgf.fwd_model= cmdl;
imgf.elem_data= imgf.elem_data(c2f_idx,:);

show_slices(imgf,levels);
print_convert dual_partial2d05a.png;

% reconstruct dual model on coarse mesh
imgd= inv_solve(drec_mdl, vh, vi);
imgd.fwd_model= cmdl;

show_slices(imgd,levels);
print_convert dual_partial2d05b.png;

