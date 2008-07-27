% Reconstruct some images $Id$
algs = {'GREIT_Sheffield_backproj', ...
        'GREIT_NOSER_ndiff', ...
        'GREIT_NOSER_diff', ...
       };

imb.calc_colours.ref_level = 0; % select colour output
imb.calc_colours.greylev   = 0.01; % black backgnd
imb.calc_colours.backgnd   = [.1,.5,.1]; %grey

idx= 157;

phi = linspace(0,2*pi,10);
xc=-xyzr_pt(1,idx) + xyzr_pt(4,idx) * sin(phi); xc= round(xc*15.5 + 16.5);
yc=-xyzr_pt(2,idx) + xyzr_pt(4,idx) * cos(phi); yc= round(yc*15.5 + 16.5);
ind= sub2ind([32,32],yc,xc);

for i= 1:length(algs)
   [img,map] = feval(algs{i}, vh, vi(:,idx) );

   imc= calc_colours(img, imb);
   imc(~map) = 1; % background
   imc(ind) = 1;

   imwrite(imc,colormap, sprintf('simulation_test03_%d.png',i),'png')
end
