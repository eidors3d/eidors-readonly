% Main code to perform analysis of system performance, as documented
%  in  the paper:  
% eidors3d.sf.net/data_contrib/evaluation-phantom/phantom_evaluation.shtml

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

[snr,accuracy,sym_err]=data_analysis(0);

RERi = recip_err;

% This line takes longer, so comment it to speed up
[drift] = drift_analysis

[DET, greit_para] = image_analysis;

%%  row 1
figure;
axes('position',[0.08 0.8 0.4 0.18])
plot(snr(:,1),'-^b','LineWidth',1,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',2);
ylabel('SNR (dB)');
xlabel('Channel number');
axis( [ 0 215 30 70 ] );

axes('position',[0.58 0.8 0.4 0.18])
plot(accuracy(:,1),'-^b','LineWidth',1,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',2);
ylabel('Accuracy');
xlabel('Channel number');
axis( [ 0 215 88 100 ] );


%   ------  row 2 -----------
axes('position',[0.09 0.55 0.4 0.18])
%ax = [1:18]*0.2;
plot( drift,'-^b','LineWidth',1,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',2);
xlabel('Averaging time (\tau)');
ylabel('Allan std');
axis( [ 0 20 0 8000 ] );

axes('position',[0.58 0.55 0.4 0.18])
plot(RERi,'-^b','LineWidth',1,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',2);
xlabel('Channel number');
ylabel('Reciprocity accuracy');
axis( [ 0 215 78 100 ] );
%  --------   row 3 -----------------

description = {'Detectability','Distinguishability'};

for i = 1:2
    if i ==1
        DET_R = DET(i,1:7);
        DET_R = fliplr(DET_R);
        DET_L = DET(i,7:13);
        ax = [0,0.14,0.29,0.43,0.57,0.71,0.86];
        axes('position',[0.08 0.3 0.4 0.18])
        plot(ax, DET_R','--rs','LineWidth',1,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor','r',...
            'MarkerSize',2);
        hold on
        plot(ax,DET_L','-bs','LineWidth',1,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor','b',...
            'MarkerSize',2);
        legend('to + x axis','to - x axis');
    else
        ax = [2,5,8,11]/14;
        axes('position',[0.58 0.3 0.4 0.18])
        plot(ax,DET(i,1:4),'-^b','LineWidth',1,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor','w',...
            'MarkerSize',2);
    end
axis( [ 0 1 22 60 ] );
xlabel('Radius fraction');
ylabel(description(i));
end

%
ax = [0,0.14,0.29,0.43,0.57,0.71,0.86];
greit_para(4,:) = greit_para(5,:); p_names={'AR','PE','RES','RNG'};
for i=1:4; %subplot(1,4,i);
    switch(i)
        case 1; axes('position',[0.07 0.04 0.18 0.18]);
        case 2; axes('position',[0.31 0.04 0.18 0.18]);
        case 3; axes('position',[0.56 0.04 0.18 0.18]);
        case 4; axes('position',[0.81 0.04 0.18 0.18]);
    end
    
    AR_R = greit_para(i,1:7);
    AR_R = fliplr(AR_R)
    AR_L = greit_para(i,7:13)
    plot(ax, AR_R','-bs','LineWidth',1,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerSize',2);
    hold on
    plot(ax,AR_L','--rs','LineWidth',1,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r',...
        'MarkerSize',2);
    ylabel(p_names{i});
    % xlabel('Radius');
    switch(i)
        case 1; axis( [ 0 1 0.7 2.5 ] );
        case 2; axis( [ 0 1 -0.05 0.5 ] );
        case 3; axis( [ 0 1 0.28 0.55 ] );
        case 4; axis( [ 0 1 2.5 4 ] );
    end
end
