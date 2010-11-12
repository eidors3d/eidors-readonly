function meas_icov = meas_icov_rm_elecs( imdl, elec_list )
% MEAS_ICOV_RM_ELECS: remove electrodes from consideration
% meas_icov = meas_icov_rm_elecs( inv_model, elec_list )
% 
% Reference Accounting for erroneous electrode data in EIT
% A. Adler Physiological Measurement, 25(1):227-238, 2004. 
%
% elec_list = numbers of elecs to remove

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

% TODO: put elec_list into a parameter
meas_icov = [];
for stim = imdl.fwd_model.stimulation(:)'
   mp = stim.meas_pattern;
   sp = stim.stim_pattern;
   icovi = ones(size(mp,1),1);
   if any(sp(elec_list) ~= 0)
      icovi = 0*icovi;
   else
      icovi = ~any( mp(:,elec_list) ~= 0, 2);
   end

   meas_icov = [meas_icov; icovi];
end
