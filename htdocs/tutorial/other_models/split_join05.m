fmdlu= join_models(fmdl1, fmdl2);
clf; subplot(121); show_fem(fmdlu); axis off; view(0,65); xlim([-1.3,1.3]);

print_convert('split_join05a.png',struct('pagesize',[12,5]));
