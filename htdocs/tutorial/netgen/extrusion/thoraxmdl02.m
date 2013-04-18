shape = { 1,                      % height
          {thorax, rlung, llung}, % contours
          [4,50],                 % perform smoothing with 50 points
          0.04};                  % small maxh (fine mesh)

elec_pos = [ 16,                  % number of elecs per plane
             1,                   % equidistant spacing
             0.5]';               % a single z-plane
         
elec_shape = [0.05,               % radius
              0,                  % circular electrode
              0.01 ]';             % maxh (electrode refinement) 

fmdl = ng_mk_extruded_model(shape, elec_pos, elec_shape);
% this similar model is also available as:
% fmdl = mk_library_model('adult_male_16el_lungs');

[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl.stimulation = stim;

img=mk_image(fmdl,1);
img.elem_data(fmdl.mat_idx{2})= 0.3; % rlung
img.elem_data(fmdl.mat_idx{3})= 0.3; % llung

clf; show_fem(img); view(0,70);
print_convert thoraxmdl02a.jpg '-density 100'
