% simulate homogeneous $Id: simulate_move2_03.m,v 1.1 2008-03-29 15:39:56 aadler Exp $

% stimulation pattern: adjacent
stim_pat= mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},1);

smdl.stimulation= stim_pat;
himg= eidors_obj('image','','fwd_model',smdl,...
                 'elem_data',ones(size(smdl.elems,1),1) );

vh= fwd_solve(himg);
