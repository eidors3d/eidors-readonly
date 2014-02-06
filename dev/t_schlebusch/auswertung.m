load('2014-02-06-sweep.mat')
%%
trainGI = data((data(:,2)==3),:)

% Show that cubic interpolation for correlation between volume and radius
% is correct (follows from V=4/3 pi r^3).
% trainGIn = trainGI;
% trainGIn(:,5) = trainGIn(:,5)-min(trainGIn(:,5));
% trainGIn(:,5) = trainGIn(:,5)./max(trainGIn(:,5));
% plot(trainGIn(:,1),trainGIn(:,5));
% hold on
% vol = (trainGIn(:,1).^3);
% vol = vol-min(vol);
% vol = vol./max(vol);
% plot(trainGIn(:,1),vol,'r');

[p_gi,ErrorEst] = polyfit(trainGI(:,5),trainGI(:,1),3);
gi_fit = polyval(p_gi,trainGI(:,5));

%%
%Plots for volume accuracy

% select simulations at v=0.5, 2, 3
tmp05 = data((data(:,2)==0.5),:);
tmp20 = data((data(:,2)==2),:);
tmp30 = data((data(:,2)==3),:);

figure
subplot(1,3,1)
hold on
plot(tmp05(:,1),tmp05(:,1),'k--');
plot(tmp05(:,1), tmp05(:,3));
plot(tmp05(:,1),polyval(p_gi,tmp05(:,5)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ v=0.5');

subplot(1,3,2)
hold on
plot(tmp20(:,1),tmp20(:,1),'k--');
plot(tmp20(:,1), tmp20(:,3));
plot(tmp20(:,1),polyval(p_gi,tmp20(:,5)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ v=2');

subplot(1,3,3)
hold on
plot(tmp30(:,1),tmp30(:,1),'k--');
plot(tmp30(:,1), tmp30(:,3));
plot(tmp30(:,1),polyval(p_gi,tmp30(:,5)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ v=3');

%%
% Plots for urine conductivity influence
tmpr01 = data((data(:,1)==0.01),:);
tmpr03 = data((data(:,1)==0.03),:);
tmpr05 = data((data(:,1)==0.05),:);

figure
subplot(1,3,1)
hold on
plot(tmpr01(:,2),tmpr01(:,1),'k--');
plot(tmpr01(:,2), tmpr01(:,3));
plot(tmpr01(:,2),polyval(p_gi,tmpr01(:,5)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ r=0.01m');
xlabel('contrast v [a.u.]');
ylabel('estimated radius r [m]');

subplot(1,3,2)
hold on
plot(tmpr03(:,2),tmpr03(:,1),'k--');
plot(tmpr03(:,2), tmpr03(:,3));
plot(tmpr03(:,2),polyval(p_gi,tmpr03(:,5)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ r=0.03m');
xlabel('contrast v [a.u.]');
ylabel('estimated radius r [m]');

subplot(1,3,3)
hold on
plot(tmpr05(:,2),tmpr05(:,1),'k--');
plot(tmpr05(:,2), tmpr05(:,3));
plot(tmpr05(:,2),polyval(p_gi,tmpr05(:,5)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ r=0.05m');
xlabel('contrast v [a.u.]');
ylabel('estimated radius r [m]');