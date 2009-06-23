% simulate targets $Id$

fmdl.stimulation= mk_stim_patterns(length(elec_nodes), 1, '{ad}','{ad}', {}, 1);

img= eidors_obj('image','fmdl','fwd_model',fmdl, ...
                'elem_data', ones(size(fmdl.elems,1),1) );
vh= fwd_solve(img);

% interpolate onto mesh
xym= interp_mesh( fmdl, 3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
% non-conductive target
ff  = (x_xym>-3) & (x_xym<-2) & (y_xym<-4) & (y_xym>-7);
img.elem_data= img.elem_data - 0.1*mean(ff,3);
% conductive target
ff  = (x_xym> 2) & (x_xym< 4) & (y_xym<-5) & (y_xym>-7);
img.elem_data= img.elem_data + 0.1*mean(ff,3);

% inhomogeneous image
vi= fwd_solve(img);

show_fem(img); axis image;
print -dpng -r125 square_mesh03a.png
