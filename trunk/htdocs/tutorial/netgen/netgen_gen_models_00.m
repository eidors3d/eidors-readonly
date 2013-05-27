netgen_gen_models;
crop_model([], inline('z<-2 | abs(x)>2.5 | abs(y)>2.5','x','y','z'));
axis([-2.5,2.5,-2.5,2.5,-2,0]);
print_convert netgen_gen_models_15a.png '-density 75';