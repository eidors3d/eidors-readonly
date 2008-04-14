% Make figures for EIDORS 3.3 paper
clf
imdm=mk_common_model('e2d4c',16);
smdl= imdm.fwd_model; % simulation model

[vh,vi,xyr_pt]= simulate_2d_movement( 20, smdl,[.75,.05]);
keep= 1:3;
vi= vi(:,keep);
xyr_pt= xyr_pt(:,keep);

% Show model and simulated targets
show_fem(smdl);
theta= linspace(0,2*pi,50); xr= cos(theta); yr= sin(theta);
hold on;
for i=1:length(xyr_pt)
    hh= plot(xyr_pt(3,i)*xr+ xyr_pt(1,i), ...
             xyr_pt(3,i)*yr+ xyr_pt(2,i));
    set(hh,'LineWidth',3,'Color',[0,0,1]);
    text(xyr_pt(1,i),xyr_pt(2,i),sprintf('%d',i), ...
        'HorizontalAlignment','center','FontSize',8, ...
        'Color',[0,0,1],'FontWeight','bold');
end
hold off;


axis image
print -dpng -r100 fig1a.png

axis image
axis([0.5,1.1,-0.5,0.1]);
print -dpng -r100 fig1b.png


imdl= mk_common_model('c2c2',16');
img= inv_solve(imdl,vh,vi);

clf; subplot(121)
img2=img;

img2.elem_data= img.elem_data(:,1);
show_fem(img2)
axis image
print -dpng -r125 fig2a.png

img2.elem_data= img.elem_data(:,2);
show_fem(img2)
axis image
print -dpng -r125 fig2b.png

img2.elem_data= img.elem_data(:,3);
show_fem(img2)
axis image
print -dpng -r125 fig2c.png

!find -name 'fig*.png' -exec convert  -trim '{}' PNG8:'{}' ';'

