xy = [0 0;  1 0; 1 1; 0 1];
mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1], [4]});
show_fem(mdl,[0 1 0]);