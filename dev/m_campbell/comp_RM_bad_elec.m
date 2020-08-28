function imdl_comp = comp_RM_bad_elec(imdl, rm_elecs, type)
% -------------------------------------------------------------------------
% Description:
%   imdl_comp = comp_RM_bad_elec(imdl, rm_elecs, type)
%
%   Modify reconstruction matrix of imdl to compensate for noisy or
%   disconnected electrodes or measurements using the method from Adler
%   2004 and Mamatjan 2013.
% -------------------------------------------------------------------------
% Parameters:
%   imdl:
%       EIDORS inverse model structure with keep intermediate results.
%   rm_elecs:
%       Electrodes to zero out from the reconstruction matrix. Can also be
%       logical matrix of measurements to zero (rather than zeroing all
%       meas from elec, just zero that bad ones).
%   type:
%       Denotes whether rm_elecs contains electrode numbers to be excluded
%       or measurements to be excluded.
% -------------------------------------------------------------------------   
% Returns:
%   imdl_comp:
%       inverse model whose reconstruction matrix has been modified with a
%       series of rank one updates that compensate for the removal of the
%       removed electrodes.
% -------------------------------------------------------------------------   
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------
% VERSION:
%   1.2.0
% -------------------------------------------------------------------------
% (C) 2019-2020 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------

imdl_comp = imdl;

if nargin == 2
    type = 'elec';
end % end if

% Find measurements from bad electrodes
if strcmp(type, 'meas')
    rm_meas = rm_elecs;
    if islogical(rm_meas)
        ee  = find(rm_meas);
    else
%         1. rmMeas  = worst(scores >= RMTHRESHOLD);
%         2. rmMeas2 = mmScores >= RMTHRESHOLD;
% order of updates affects results. Sort measurements for consistency
% between 1.  measurements selected for rejection from array of
% measurements sorted by badness and 2. from array of all measurements
% greater than threshold.
        ee = sort(rm_meas); 
    end % end if
else
    kk = meas_icov_rm_elecs(imdl_comp, rm_elecs);
    ee = sort(find(diag(kk)~=1)); % bad channels of meas_sel
end % end if

% do work on reconstruction matrix
imdl_comp.solve_use_matrix.RM_orig = imdl_comp.solve_use_matrix.RM;
PJt     = imdl_comp.solve_use_matrix.PJt;
X       = imdl_comp.solve_use_matrix.X;
X_star  = X;

for i= 1:length(ee)
    channel             = ee(i);
    Xi                  = X_star(:, channel);
    X_star              = X_star- Xi* Xi' * inv(X_star(channel, channel)); % update on X_star each time. go back and explicitly zero bad meas
    X_star(channel, :)  = 0;
    X_star(:, channel)  = 0;
end % end for

imdl_comp.solve_use_matrix.X_star   = X_star;
imdl_comp.solve_use_matrix.RM       = PJt* X_star;

end % end function