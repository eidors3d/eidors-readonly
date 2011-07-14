
% create model $Id$
%% 1
% Model parameters
n_elec= 16;
n_nodes= 2000;

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

%% 2.
% Create an object within the model
trg_ctr= [.7,.1];
trg_rad= 0.1;
th=linspace(0,2*pi,20)';th(end)=[];
trg_refine_nodes= [refine_nodes; ...
             trg_rad*sin(th)+trg_ctr(1), ...
             trg_rad*cos(th)+trg_ctr(2) ];

tmdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, elec_nodes, ...
                       trg_refine_nodes, z_contact);

% find elements in size target
mdl_pts = interp_mesh( tmdl );
ctr_pts = mdl_pts - ones(size(mdl_pts,1),1)*trg_ctr;
in_trg  =(sum(ctr_pts.^2,2) < trg_rad^2);


timg= eidors_obj('image','','fwd_model',tmdl,...
                 'elem_data',1 + in_trg*.5);

% % Show output - full size
% subplot(121); show_fem( smdl ); axis equal; axis([-1.1 1.1 -1.1 1.1]);  
% subplot(122); show_fem( timg ); axis equal; axis([-1.1 1.1 -1.1 1.1]);
% print_convert simulate_move2_02c.png '-density 125'
% 
% % Show output - full size
% subplot(121); show_fem( smdl ); axis([.5,1.05,-.1,.3]); axis on; axis equal
% subplot(122); show_fem( timg ); axis([.5,1.05,-.1,.3]); axis on; axis equal
% print_convert simulate_move2_02c.png '-density 125'

%% 3 simulate homogeneous $Id$

% stimulation pattern: adjacent
stim_pat= mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},1);

smdl.stimulation= stim_pat;
himg= eidors_obj('image','','fwd_model',smdl,...
                 'elem_data',ones(size(smdl.elems,1),1) );

vh= fwd_solve(himg);

%% 4 moving target $Id$

% Create a moving object within the model
trg_rad= 0.1;
radius= 0.75;
n_sims= 20;
contrast= 0.1;

clear vi;
for i= 1:n_sims;
   thc= 2*pi*(i-1)/n_sims;
   trg_ctr= radius*[cos(thc),sin(thc)];

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

   clf; show_fem(timg); axis equal
   print_convert(sprintf('simulate_move2_04c%02d.png',i),'-density 50');
   vi(i)= fwd_solve(timg);
end

%% 5 Create animated graphics $Id$

% Trim images
!find -name 'simulate_move2_04c*.png' -exec convert  -trim '{}' PNG8:'{}' ';'

% Convert to animated Gif
!convert -delay 50 simulate_move2_04c*.png -loop 0 simulate_move2_05c.gif

%% 6 Reconstruct images $Id$

imdl= mk_common_model('c2c2',16);
img= inv_solve(imdl, vh, vi);
subplot(121); show_fem(img); axis image
subplot(122); show_slices(img)

print_convert simulate_move2_06c.png '-density 125'
