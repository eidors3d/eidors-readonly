% Sample Data $Id$

if ~exist('simulate_2d_movement_data.mat')
    [vh,vi,xyr_pt]=simulate_2d_movement;
    save simulate_2d_movement_data vi vh xyr_pt
else
    load simulate_2d_movement_data
end

% Temporal solver works best at large Noise
nsr= 4.0; %Noise to Signal Ratio

signal = mean( abs(mean(vi,2) - vh) ); % remember to define signal
% Only add noise to vi. This is reasonable, since we have
% lots of data of vh to average to reduce noise
vi_n= vi + nsr*signal*randn( size(vi) );
