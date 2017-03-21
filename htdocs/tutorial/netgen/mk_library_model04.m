fmdl= mk_library_model('adult_male_grychtol2016a_2x16')
subplot(211);
view(-20,24); show_fem_enhanced(fmdl);
print_convert mk_library_model04a.png
view(0,0);    show_fem_enhanced(fmdl);
print_convert mk_library_model04b.png
