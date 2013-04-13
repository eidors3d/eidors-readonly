img_v = img;
% Stimulate between elecs 16 and 5 to get more interesting pattern
img_v.fwd_model.stimulation(1).stim_pattern = sparse([16;5],1,[1,-1],16,1);
img_v.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img_v);

img_v = rmfield(img, 'elem_data');
img_v.node_data = vh.volt(:,1);
img_v.calc_colours.npoints = 128;

PLANE= [inf,inf,0.35]; % show voltages on this slice

subplot(221);
show_slices(img_v,PLANE); axis off; axis equal
print_convert thoraxmdl03a.jpg
% 
