function example_images( basename, alg_names )
% Example GREIT images $Id$
[s_algs,f_algs]= regexp(alg_names,'(\S+)')

load testdata.mat

imb.calc_colours.ref_level = 0; % select colour output
imb.calc_colours.greylev   = 0.01; % black backgnd
imb.calc_colours.backgnd   = [.5,.5,.5]; %grey

[x,y]= meshgrid(linspace(-1,1,32),linspace(-1,1,32));
map = x.^2 + y.^2 < 1.1;

for i= 1:length(s_algs)-1 
   fname = alg_names( s_algs(i+1) : f_algs(i+1));
   load(fname);

   for k= 1:length(test_v)
      dd= test_v(k);
    
      if normalize_flag
         dv = ( dd.vi - dd.vh ) ./ dd.vh;
      else
         dv = ( dd.vi - dd.vh );
      end
      ds = RM*dv;

      img= reshape(ds, 32,32);
   
      imc= calc_colours(img, imb);
      imc(~map) = 1; % background

      imwrite(imc,colormap, sprintf([ ...
        basename, '_%d%d.png'],i,k),'png')
    end
end
