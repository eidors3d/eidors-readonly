% simulate targets $Id$

fmdl.stimulation= mk_stim_patterns(n_elec, 1, [0,1],[0,1], {}, 1);

img= eidors_obj('image','fmdl','fwd_model',fmdl, ...
                'elem_data', ones(size(fmdl.elems,1),1) );
vh= fwd_solve(img);

% interpolate onto mesh
xym= interp_mesh( fmdl, 3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff  = (x_xym>-3) & (x_xym<-2) & (y_xym<-2) & (y_xym>-3);
ff = mean(ff,3);
img.elem_data= img.elem_data + ff;

% inhomogeneous image
vi= fwd_solve(img);

show_fem(img); axis image;
print_convert dm_geophys03a.png '-density 125'
