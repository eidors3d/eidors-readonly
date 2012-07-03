% moving target $Id$

% Create a moving object within the model
trg_rad= 0.05; trs= trg_rad*sin(th); trc= trg_rad*cos(th);
radius= 0.9;
n_sims= 2 ;
contrast= 0.1;

vi=[];
for i= 1:n_sims; 
   f_frac = (i-1)/(n_sims-1);
   cv= 2*pi*f_frac * 73;
   trg_ctr= f_frac*radius*[cos(cv),sin(cv)];

   trg_refine_nodes= [refine_nodes; [ trs+trg_ctr(1), trc+trg_ctr(2) ]];

   tmdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, elec_nodes, trg_refine_nodes, z_contact);
   tmdl.stimulation = stim_pat;

   % find elements in size target
   mdl_pts = interp_mesh( tmdl );
   ctr_pts = mdl_pts - ones(size(mdl_pts,1),1)*trg_ctr;
   in_trg  =(sum(ctr_pts.^2,2) < trg_rad^2);

   % Create target image object
   timg= eidors_obj('image','','fwd_model',tmdl, 'elem_data',1 + in_trg*contrast);
   vi_t= fwd_solve(timg);
   vi = [vi, vi_t.meas];
   if     i==1;    show_fem(timg); print -dpng -r60 simulate_2d_test03a.png
   elseif i==150;  show_fem(timg); print -dpng -r60 simulate_2d_test03b.png
   end
end
