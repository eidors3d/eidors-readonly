% Lung images
% $Id$

% 2D Model
imdl= mk_common_model('c2t3',16);


% most EIT systems image best with normalized difference
imdl.fwd_model.normalize_measurements= 1;
imdl.RtR_prior= @prior_gaussian_HPF;

% electrodes start on back (dorsal), then do this
imdl.fwd_model.electrode([9:16,1:8])=  ...
   imdl.fwd_model.electrode;


subplot(221);
show_fem(imdl.fwd_model);

axis equal
print_convert tutorial410a.png;
