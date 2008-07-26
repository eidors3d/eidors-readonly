% 2D image $Id$

img= eidors_obj('image','2D rectangle', ...
      'elem_data', ones(size(mdl.elems,1),1) * conduc, ...
      'fwd_model', mdl); 

fsol= fwd_solve(img)
%show_fem(mdl); axis('equal'); set(gca,'Ylim',[-.5,ww-.5]);
%print -r100 -dpng tutorial022b.png
