fmdlu= join_models(fmdl1, fmdl2);

clf; subplot(131); show_fem(fmdlu, [0,1,1]); axis off;

print_convert('split_join03a.png',struct('pagesize',[15,5]));
