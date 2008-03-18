% Dual Partial $Id: dual_partial2d04.m,v 1.1 2008-03-18 18:59:47 aadler Exp $

clf; levels= [inf,inf,0,1,1];

% reconstruct fine model
imgf= inv_solve(frec_mdl, vh, vi);

show_slices(imgf,levels);
print -r125 -dpng dual_partial2d04a.png;

% reconstruct dual model
imgd= inv_solve(drec_mdl, vh, vi);

show_slices(imgd,levels);
print -r125 -dpng dual_partial2d04b.png;

