% Demo algorithms $Id$

% SELECT ALGORITHMS
algs = {'GREIT_Sheffield_backproj', ...
        'GREIT_NOSER_ndiff', ...
        'GREIT_NOSER_diff', ...
       };

% SELECT COLOUR OUTPUT
imb.calc_colours.ref_level = 0;
imb.calc_colours.greylev   = 0.01; % black backgnd
imb.calc_colours.backgnd   = [.5,.5,.5]; %grey


for i= 1:length(algs)
   for k= 1:4
      [img,map] = feval(algs{i}, v(k).vh, v(k).vi );

      imc= calc_colours(img, imb);
      imc(~map) = 1; % background

      imwrite(imc,colormap, ...
              sprintf('demo_algs04_%d%d.png',i,k),'png')
    end
end
