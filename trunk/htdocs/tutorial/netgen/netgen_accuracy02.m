maxsz = [8,4,3,2,1.5,1.3,1.2,1.1,1,.9,0.8,0.7,0.6,0.5,0.45,0.4,0.35,0.3,0.28,0.26,.25,.24,.23,.22,.21,.20,.19,.18,.17,.16];
%maxsz = [8,4,3,2,1,0.8,0.7];
for i=1:length(maxsz)
   fmdl = ng_mk_cyl_models([2*h,w/2,maxsz(i)],[16,h],elec_sz); 
   fmdl.stimulation = stim;
   vh = fwd_solve(mk_image(fmdl,1));
   vh.n_ne = [size(fmdl.nodes,1), size(fmdl.elems,1)];
   vv(i) = vh;
end
