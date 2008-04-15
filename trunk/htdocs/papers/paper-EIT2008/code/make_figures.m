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
imdl.RtR_prior= @gaussian_HPF_prior;
iml.hyperparameter.value= 0.003;
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

% MOVING BALL #2
% Model parameters
n_elec= 16;
n_nodes= 5000;

% Create electrodes
refine_level=4; %electrode refinement level
elec_width= .1;
z_contact = 0.01;
% elect positions
th=linspace(0,2*pi,n_elec(1)+1)';th(end)=[];
elec_posn= [sin(th),cos(th)];
[elec_nodes, refine_nodes] = dm_mk_elec_nodes( elec_posn, ...
       elec_width, refine_level);

% Define circular medium
fd=inline('sum(p.^2,2)-1','p');
bbox = [-1,-1;1,1];
smdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, ...
                          elec_nodes, refine_nodes, z_contact);

% stimulation pattern: adjacent
stim_pat= mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},1);

smdl.stimulation= stim_pat;
himg= eidors_obj('image','','fwd_model',smdl,...
                 'elem_data',ones(size(smdl.elems,1),1) );

vh= fwd_solve(himg);

% Create a moving object within the model
trg_rad= 0.1;
radius= 0.75;
n_sims= 20
contrast= 0.1;
if n_nodes>2000
   n_th_obj= 30;
else
   n_th_obj= 15;
end

th=linspace(0,2*pi,n_th_obj+1)';th(end)=[];
clear vi;
for i= 1:3;
   thc= 2*pi*(i-1)/n_sims;
   trg_ctr= radius*[cos(thc),-sin(thc)];

   trg_refine_nodes= [refine_nodes; ...
                trg_rad*sin(th)+trg_ctr(1), ...
                trg_rad*cos(th)+trg_ctr(2) ];

   tmdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, elec_nodes, ...
                          trg_refine_nodes, z_contact);
   tmdl.stimulation = stim_pat;

   % find elements in size target
   mdl_pts = interp_mesh( tmdl );
   ctr_pts = mdl_pts - ones(size(mdl_pts,1),1)*trg_ctr;
   in_trg  =(sum(ctr_pts.^2,2) < trg_rad^2);

   % Create target image object
   timg= eidors_obj('image','','fwd_model',tmdl,...
                    'elem_data',1 + in_trg*contrast);
   vi(i)= fwd_solve(timg);

   if i==1
      clf; show_fem(timg);
      axis image
      print -dpng -r100 fig3a.png
      
      axis image
      axis([0.5,1.1,-0.5,0.1]);
      print -dpng -r100 fig3b.png

   end
end


show_fem(inv_solve(imdl,vh,vi(1)));
axis image; print -dpng -r125 fig4a.png
show_fem(inv_solve(imdl,vh,vi(2)));
axis image; print -dpng -r125 fig4b.png
show_fem(inv_solve(imdl,vh,vi(3)));
axis image; print -dpng -r125 fig4c.png
!find -name 'fig*.png' -exec convert  -trim '{}' PNG8:'{}' ';'

