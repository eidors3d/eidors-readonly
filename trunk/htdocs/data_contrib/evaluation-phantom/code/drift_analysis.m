function [drift] = drift_analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Name:drift_analysis.m
% Comments:
% input: nframes, nmeas
% nframes is the number of frames for drift measurement
% nmeas is number of measurements for each frame (mostly 208, 416)
% NOTE: the different number of measurements will affact Allan std..
% recommend 208 for the consistency of preseting results
% output: sig is the values of Allan variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

nmeas = 208; %mnumber of easurements for each frame (mostly 208)
nframes=300; % number of frames for drift measurement
Nseq = 1;   % number of sequence that taken

% load measured data and arrange it in an array in a sequence
VH = load_meas(nmeas, nframes, Nseq);

%% calculate Allan std
tau0 =  nmeas;
[drift,sig2,osig,msig,tsig,tau]=avar(VH,tau0);

% plot the Allan std
plotd(drift)
end

function VH=load_meas(nmeas, nframes, Nseq)

% Load measurement files
NPos = {'m1.mat';'m2.mat';'m3.mat';'m4.mat';'m5.mat';'m6.mat';'m7.mat';'m8.mat'; ...
    'm9.mat'; 'm10.mat';'m11.mat';'m12.mat';};
datapath =  ('.././data/drift/2010/');
IDX = [18:2:416,2:2:16];

Tot = nframes*nmeas;
% rearrange data into one array
vh =0; VH = [];
for ii = 1 : Nseq
    vh = load([datapath, NPos{ii}]);
    v1 = [vh.Eit_Data{:}];
    v2= abs(v1(IDX,:));
    v= reshape(v2,Tot,1);
    VH = [VH; v];
end
end

function plotd(drift)
plot(drift,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',4);
ylabel(' Allan std');
xlabel('Averaging time (\tau)');
% axis( [ 0 20 0 7500 ] );
end
