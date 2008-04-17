% simulate targets $Id: square_mesh03.m,v 1.1 2008-04-17 16:12:05 aadler Exp $

fmdl.stimulation= mk_stim_patterns(n_elec, 1, '{ad}','{ad}', {}, 1);

img= eidors_obj('image','fmdl','fwd_model',fmdl, ...
                'elem_data', ones(size(fmdl.elems,1),1) );
vh= fwd_solve(img);

% interpolate onto mesh
xym= interp_mesh( fmdl, 3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff  = (x_xym>-3) & (x_xym<-2) & (y_xym<-2) & (y_xym>-5);
ff = mean(ff,3);
img.elem_data= img.elem_data + ff;

% inhomogeneous image
vi= fwd_solve(img);

show_fem(img); axis image;
print -dpng -r125 square_mesh03a.png
