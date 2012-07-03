fmdlr = ng_mk_extruded_model({0,trunk,[4,50],.1},[0,0],[0.1]);
fmdlr.nodes = fmdlr.nodes*diag([-1,-1]);

show_fem(fmdlr); view(0,90);
print_convert pig_body07.jpg
