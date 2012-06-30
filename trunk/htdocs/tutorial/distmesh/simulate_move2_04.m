% moving target $Id$

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
   timg= mk_image(tmdl,1 + in_trg*contrast);

   clf; show_fem(timg); axis equal
   print_convert(sprintf('simulate_move2_04a%02d.png',i),'-density 50');
   vi(i)= fwd_solve(timg);
end
