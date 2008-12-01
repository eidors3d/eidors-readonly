function print_thorax_model( printfile, model )
   sp = find (model==' '); 
   load(model(sp+1:end));
   show_fem(fmdl);
   view(35,10)
   print('-dpng','-r150',printfile);
