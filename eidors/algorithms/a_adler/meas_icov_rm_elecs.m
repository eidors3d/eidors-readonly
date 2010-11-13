function meas_icov = meas_icov_rm_elecs( imdl, elec_list )
% MEAS_ICOV_RM_ELECS: remove electrodes from consideration
% meas_icov = meas_icov_rm_elecs( inv_model, elec_list )
%
% elec_list = numbers of elecs to remove
%   or as imdl.meas_icov_rm_elecs.elec_list
% 
% Reference Accounting for erroneous electrode data in EIT
% A. Adler Physiological Measurement, 25(1):227-238, 2004. 

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin< 2
   elec_list = imdl.meas_icov_rm_elecs.elec_list;
end

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

n = length(meas_icov);
meas_icov = spdiags( meas_icov, 0, n,n );
