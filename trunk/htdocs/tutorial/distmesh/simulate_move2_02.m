% create model $Id$

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

% Show output - full size
subplot(121); show_fem( smdl ); axis on;
subplot(122); show_fem( timg ); axis on;
print -dpng -r125 simulate_move2_02a.png

% Show output - full size
subplot(121); show_fem( smdl ); axis([.5,1.05,-.1,.3]); axis on;
subplot(122); show_fem( timg ); axis([.5,1.05,-.1,.3]); axis on;
print -dpng -r125 simulate_move2_02b.png
