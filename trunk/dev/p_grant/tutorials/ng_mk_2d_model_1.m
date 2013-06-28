xy = [0 0;  1 0; 1 1; 0 1];
mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [0.5 1; 0.5 0; 0 0.5]);
show_fem(mdl,[0 1 0]);