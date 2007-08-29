function transfimp = calc_transferimpedance( img)
% CALC_TRANSFERIMPEDANCE: Calculates transfer impedance
% 
%   transfimp = calc_transferimpedance( fwd_model, img)
%
% fwd_model is a fwd_model structure
% img       is an image structure
% (C) 2006 Bill Lionheart. License: GPL version 2 or version 3
% $Id: calc_transferimpedance.m,v 1.1 2007-08-29 09:07:17 aadler Exp $

% create new stim patterns
% stimulate with one ref electrode and then each in turn
% make an equiv set of measurements

n_elecs= length( img.fwd_model.electrode );

%[stim_pat, meas_pat]= trigonometric( n_elecs );
%[stim_pat, meas_pat]= electrode_wise( n_elecs );
%[stim_pat, meas_pat]= monopolar( n_elecs );
 [stim_pat, meas_pat]= monopolar_even( n_elecs );
img.fwd_model.stimulation = stim_pat;

lambda= 1e-3;
imeas_pat= pinv(meas_pat);

data = fwd_solve(img);

sz= length(img.fwd_model.stimulation);
transfimp = reshape( data.meas, sz, sz);
transfimp = imeas_pat * transfimp * imeas_pat';

function [stim_pat, meas_pat] = trigonometric( n_elecs )
    stim_pat = struct;
    idx= linspace(0,2*pi,n_elecs+1)'; idx(end)= [];
    omega= idx*[1:n_elecs/2];
    meas_pat= [cos(omega), sin(omega) ]';
    for i=1:n_elecs
        stim_pat(i).stimulation='mA';
        stim_pat(i).stim_pattern= meas_pat(i,:)';
        stim_pat(i).meas_pattern= meas_pat;
    end

function [stim_pat, meas_pat] = electrode_wise( n_elecs)
    stim_pat = struct;
    meas_pat= [-ones(n_elecs-1,1), speye(n_elecs-1)];
    for i=2:n_elecs
        stim_pat(i-1).stimulation='mA';
        stim_pat(i-1).stim_pattern= sparse([1,i],1,[-1,1],n_elecs,1);
        stim_pat(i-1).meas_pattern= meas_pat;
    end

function [stim_pat, meas_pat] = monopolar( n_elecs)
    stim_pat = struct;
    meas_pat= speye(n_elecs);
    for i=1:n_elecs
        stim_pat(i).stimulation='mA';
        stim_pat(i).stim_pattern= sparse(i,1,1,n_elecs,1);
        stim_pat(i).meas_pattern= meas_pat;
    end

function [stim_pat, meas_pat] = monopolar_even( n_elecs)
    stim_pat = struct;
    meas_pat= eye(n_elecs) - ones(n_elecs)/n_elecs;
    for i=1:n_elecs
        stim_pat(i).stimulation='mA';
        stim_pat(i).stim_pattern= meas_pat(i,:)';
        stim_pat(i).meas_pattern= meas_pat;
    end
