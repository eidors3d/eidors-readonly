function example_images( basename )
% Example GREIT images $Id$


algs= get_list_of_algs;
load testdata.mat

imb.calc_colours.ref_level = 0; % select colour output
imb.calc_colours.greylev   = 0.01; % black backgnd
imb.calc_colours.backgnd   = [.5,.5,.5]; %grey

for i= 1:length(algs)
   for k= 1:4
      [img,map] = feval(algs{i}, test_v(k).vh, test_v(k).vi );

      imc= calc_colours(img, imb);
      imc(~map) = 1; % background

      imwrite(imc,colormap, sprintf([ ...
        basename, '_%d%d.png'],i,k),'png')
    end
end
