% Create model for Absolute Imaging
% Based on example from Stephen Murphy
% $Id$

load tutorial151_model.mat

% Create Simulation Model
sim_img= eidors_obj('image', 'Simulation Image');
sim_img.fwd_model= simmdl;

backgnd= .02; 
sim_img.elem_data= backgnd*ones(size(simmdl.elems,1),1);
sim_img.elem_data(target1)= backgnd*2;
sim_img.elem_data(target2)= backgnd/4;

clf;
axes('position',[.1,.1,.65,.8]);
   show_fem(sim_img); view(-33,20);
axes('position',[.8,.1,.15,.8]);
   show_slices(sim_img,5);

print_convert tutorial151a.png '-density 75'
