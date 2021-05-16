function imdl_comp = comp_RM_bad_elec(imdl, rm_data, type)
% -------------------------------------------------------------------------
% Description:
%   imdl_comp = comp_RM_bad_elec(imdl, rm_data, type)
%
%   Modify reconstruction matrix of imdl to compensate for noisy or
%   disconnected electrodes or measurements using the method from Adler
%   2004 and Mamatjan 2013.
% -------------------------------------------------------------------------
% Parameters:
%   imdl (EIDORS inverse model object):
%       EIDORS inverse model structure with keep intermediate results.
%   rm_data (array):
%       Electrodes to zero out from the reconstruction matrix. Can also be
%       logical matrix of measurements to zero (rather than zeroing all
%       meas from elec, just zero that bad ones).
%   type (string):
%       Denotes whether rm_elecs contains electrode numbers to be excluded
%       or measurements to be excluded.
% -------------------------------------------------------------------------   
% Returns:
%   imdl_comp (EIDORS inverse model object):
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

% Find measurements from bad electrodes
if strcmp(type, 'meas')
    rm_meas = rm_data;
    if islogical(rm_meas)
        ee = find(rm_meas);
    else
        ee = sort(rm_meas); 
    end % end if
else
    ee = electrodes_used_meas_or_stim(rm_data, imdl);
    ee = sort(ee);
end % end if

% do work on reconstruction matrix
imdl_comp.solve_use_matrix.RM_orig = imdl_comp.solve_use_matrix.RM;
PJt = imdl_comp.solve_use_matrix.PJt;
X = imdl_comp.solve_use_matrix.X;
X_star = X;

for i = 1: length(ee)
    channel = ee(i);
    Xi = X_star(:, channel);
    % update on X_star; go back and explicitly zero bad meas
    X_star = X_star - Xi * Xi' * inv(X_star(channel, channel));
    X_star(channel, :) = 0;
    X_star(:, channel) = 0;
end % end for

imdl_comp.solve_use_matrix.X_star = X_star;
imdl_comp.solve_use_matrix.RM = PJt * X_star;

end % end function


function wasUsed = electrodes_used_meas_or_stim(rm_elecs, imdl)
    nElecs = length(imdl.fwd_model.electrode);
    nMeas = sum(imdl.fwd_model.meas_select);    
    wasUsed = zeros(nMeas, 1);
    measPerPair = nMeas / nElecs;
    
    count = 1;
    for i = 1: nElecs
        [r, c] = find(imdl.fwd_model.stimulation(i).meas_pattern);
        for j = 1:max(r)
            temp = ismember(rm_elecs, c(r == j)');
            wasUsed(count) = sum(temp) > 0;
            count = count + 1;
        end % end for
    end % end for
    for i = 1: length(rm_elecs)
        thisElec = rm_elecs(i);
        stop  = thisElec * measPerPair;
        start = stop - measPerPair + 1;
        wasUsed(start : stop) = 1;
    end
    wasUsed = find(wasUsed);
end