% Reconstruct some images $Id$
load sim_radmove_homog.mat

imb.calc_colours.ref_level = 0; % select colour output
imb.calc_colours.greylev   = 0.01; % black backgnd
bkgnd= [.1,.5,.1]; imb.calc_colours.backgnd   = bkgnd;

idx= 157; dir='simulation_test_imgs';

phi = linspace(0,2*pi,10);
xc=-xyzr_pt(1,idx) + xyzr_pt(4,idx) * sin(phi); xc= round(xc*15.5 + 16.5);
yc=-xyzr_pt(2,idx) + xyzr_pt(4,idx) * cos(phi); yc= round(yc*15.5 + 16.5);
ind= sub2ind([32,32],yc,xc);

% Use this Map for reconstruction shape
[x,y]=meshgrid(linspace(-1,1,32),linspace(-1,1,32)); out = x.^2+y.^2>1.1;

algs = get_list_of_algs;
for i= 1:length(algs)
   img = feval(algs{i}, vh, vi(:,idx) );
   hmi = calc_hm_set( img, 0.5 )+1;
   qmi = calc_hm_set( img, 0.25 )+1;

   imc= calc_colours(img, imb);
   imc(out) = 1; hmi(out) = 3; qmi(out) = 3; % background
   imc(ind) = 1; hmi(ind) = 3; qmi(ind) = 3; % target

   imwrite(imc,colormap, sprintf('%s/simulation_test03_%d.png',dir,i),'png')
   clrmap = [0,0,0;1,1,1;bkgnd];
   imwrite(hmi,clrmap, sprintf('%s/simulation_test03_h%d.png',dir,i),'png')
   imwrite(qmi,clrmap, sprintf('%s/simulation_test03_q%d.png',dir,i),'png')
end
