function [DET, greit_para] = image_analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:image_analysis.m
% image_analysis include detectability, GREIT measures and distinguishability
% Comments: detectability analysis - one objects moves horizontally from the edge to
% centeral position of the tank and then to the edge
% Distinguishability analysis-two objects moving away from the centeral
% position. This algorithm calculates both in a loop: 1for both detectability;  
% 2 for distinguishability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

% initialization
description = {'Detectability vlaue','Distinguishability vlaue'};
Dist =[1,2];   % 1for both detectability;  2 for distinguishability
GREIT_para_Def =[1,0];% for 1-object detectability, not for 2-object distinguishability
Obj_radius = 2.3/14; % normalized object radius with the tank radius

for jj = Dist
    % load parameters and meas data
    [vh, vi, NumPos, XPos] = load_meas_data(jj);
    
    %% image reconstruction
    [imgr, imgav]=calc_inverse( vh, vi, NumPos);
    
    %% calc detectability or distinguishability values
    DET(jj,1:NumPos) = calc_distinguish(NumPos,imgr);
    plot_img(DET(jj,1:NumPos),description(jj),imgav,XPos)
    
    %% calculate GREIT measures
    if GREIT_para_Def(jj) % for 1-object  detectability, not for 2-object distinguishability
        figure
        greit_para=GREIT_params(imgav,XPos,NumPos,Obj_radius);
    end
end
end

function [vh, vi, NumPos, XPos] = load_meas_data(jj)
% Load measurement sequences for target movement
NPos = {'m1.mat';'m2.mat';'m3.mat';'m4.mat';'m5.mat';'m6.mat';'m7.mat';'m8.mat'; ...
    'm9.mat'; 'm10.mat';'m11.mat';'m12.mat';'m13.mat';};

switch jj
    case 1; % for detectability measure
        NumPos = 13;
        datapath =  ('.././data/detect1obj/');
        XPos = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]/14; % normalized object positions
        vh = load([datapath,'homo.mat']);
    case 2;% for distinguishability measure
        NumPos = 4;
        datapath =  ('.././data/distinguish2obj/');
        XPos = [2,5,8,11]/14; % normalized target positions in the tank
        vh = load([datapath,'obj1.mat']);
    otherwise; error('???');
end

for ii = 1 : NumPos
    vi{ii} = load([datapath,NPos{ii}]);
end
end


function [imgr, imgav]=calc_inverse(vh, vi, NumPos)
% Create Model
imdl=mk_common_model('d2c2',16);

% Set Parameters
imdl.hyperparameter.value = .15;
imdl.fwd_model.normalize_measurements = 1;

% Obtain FILTERED JACOBIAN / Protocol
% required for system A only, comment it if it is not necessary
imdl = A_filter_protocol(imdl);
vh = abs([vh.Eit_Data{:}]);
% inverse model
for ii = 1 : NumPos
    vii = abs([vi{ii}.Eit_Data{:}]);
    % Solve inverse model
    imgr(ii)= inv_solve(imdl, vh, vii);
    imgav(ii) = inv_solve(imdl, mean(vh,2), mean(vii,2));
end
end

function Det = calc_distinguish(NumPos,imgr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:calc_distinguish.m
% Comments:
% NumPos is the number of positions and measurements. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : NumPos
    % Detectability evaluation - ROI
    VOL = get_elem_volume(imgr(1).fwd_model);
    ROI = ones(size(VOL));
    avg_ed = mean(imgr(ii).elem_data,2);
    [maxval,maxpt] = max(abs(avg_ed));
    qmax = avg_ed(maxpt)/4;
    ROI( sign(qmax)*avg_ed < abs(qmax) ) = 0;
    
    VOL = VOL.*ROI;
    
    sig(ii,:) = VOL'*imgr(ii).elem_data(:,:);
    mn(ii)= mean(sig(ii,:)); st(ii)= std(sig(ii,:)); zy_ROI(ii) = mn(ii)/st(ii);
    Det(ii)=log10(abs(mn(ii)/st(ii)))*20;
end
end

function params=GREIT_params(imgav,XPos,NumPos,Obj_radius)

YPos=zeros(1,NumPos); ZPos=zeros(1,NumPos);
xyzr_pt = [XPos;YPos;ZPos;ones(1,NumPos)*Obj_radius];

imgav(1).elem_data = [imgav(:).elem_data];
params = eval_GREIT_fig_merit(imgav(1), xyzr_pt);

p_names = {'AR','PE','RES','SD','RNG'};

for i=1:5; subplot(5,1,i);
    plot(XPos,params(i,:)); ylabel(p_names{i});
    if i<5; set(gca,'XTickLabel',[]);end
end
xlabel('Radius fraction');
end

function plot_img(Det,description,imgav,XPos)
% plot detectability or distinguishability values
figure
plot(XPos,Det,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',4);

xlabel('Radius fraction');
ylabel(description);
% axis( [ -10 10 15 50 ] );

% plotting reconstructed image
imgav(1).elem_data = [imgav(:).elem_data];
% figure; show_slices(imgav(1));
end
