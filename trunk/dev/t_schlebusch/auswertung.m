load('ws-2014-02-09--big-sim-results.mat')
%%
data = cell2mat(results);
trainGI = data((data(:,2)==2),:) % select contrast=2 as reference

[p_gi,ErrorEst] = polyfit(trainGI(:,6),trainGI(:,1),3);
gi_fit = polyval(p_gi,trainGI(:,6));

%%
%Plots for volume accuracy

% select simulations at v=1, 2, 3
tmp05 = data((data(:,2)==1),:);
tmp05 = tmp05((tmp05(:,3)==1e4),:);
tmp20 = data((data(:,2)==2),:);
tmp20 = tmp20((tmp20(:,3)==1e4),:);
tmp30 = data((data(:,2)==3),:);
tmp30 = tmp30((tmp30(:,3)==1e4),:);

figure
subplot(1,3,1)
hold on
plot(tmp05(:,1),tmp05(:,1),'k--');
plot(tmp05(:,1), tmp05(:,4));
plot(tmp05(:,1),polyval(p_gi,tmp05(:,6)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ v=0.5');
ylimits = [0.2 0.5]

subplot(1,3,2)
hold on
plot(tmp20(:,1),tmp20(:,1),'k--');
plot(tmp20(:,1), tmp20(:,4));
plot(tmp20(:,1),polyval(p_gi,tmp20(:,6)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ v=2');

subplot(1,3,3)
hold on
plot(tmp30(:,1),tmp30(:,1),'k--');
plot(tmp30(:,1), tmp30(:,4));
plot(tmp30(:,1),polyval(p_gi,tmp30(:,6)),'r')
legend('ideal','parametric','GI');
title('parametric vs. GI @ v=3');

%%
% Error plots
% select simulations at v=1, 2, 3
tmp05 = data((data(:,2)==1),:);
tmp05 = tmp05((tmp05(:,3)==1e4),:);
tmp20 = data((data(:,2)==2),:);
tmp20 = tmp20((tmp20(:,3)==1e4),:);
tmp30 = data((data(:,2)==3),:);
tmp30 = tmp30((tmp30(:,3)==1e4),:);

figure
hold on
plotthis=tmp05;
plot(plotthis(:,1),plotthis(:,1)-plotthis(:,1),'k--');
plot(plotthis(:,1), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'+');
plot(plotthis(:,1),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'o')

plotthis=tmp20;
plot(plotthis(:,1), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'x');
plot(plotthis(:,1),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'s')

plotthis=tmp30;
plot(plotthis(:,1), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'d');
plot(plotthis(:,1),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'v')

legend('ideal','parametric v=1','GI v=1','parametric v=2','GI v=2','parametric v=3','GI v=3');
title('parametric vs. GI');
ylabel('relative error')
xlabel('bladder radius [m]');

%%
% Conductivity influence
% select simulations at different volumes
tmp05 = data((data(:,1)<0.02286),:);
tmp05 = tmp05((tmp05(:,3)==1e4),:);
tmp20 = data((data(:,1)>0.0390),:);
tmp20 = tmp20((tmp20(:,1)<0.0392),:);
tmp20 = tmp20((tmp20(:,3)==1e4),:);
tmp30 = data((data(:,1)>0.0507),:);
tmp30 = tmp30((tmp30(:,3)==1e4),:);

figure
hold on
plotthis=tmp05;
plot(plotthis(:,2),plotthis(:,1)-plotthis(:,1),'k--');
plot(plotthis(:,2), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'+');
plot(plotthis(:,2),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'o')

plotthis=tmp20;
plot(plotthis(:,2), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'x');
plot(plotthis(:,2),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'s')

plotthis=tmp30;
plot(plotthis(:,2), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'d');
plot(plotthis(:,2),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'v')

legend('ideal','parametric 50 ml','GI 50 ml','parametric 250 ml','GI 250 ml','parametric 550 ml','GI 550 ml');
title('parametric vs. GI');
ylabel('relative error')
xlabel('bladder contrast [a.u.]');

%%
% Noise stability
% select simulations at different volumes
tmp05 = data((data(:,1)<0.02286),:);
tmp05 = tmp05((tmp05(:,2)==2),:);
tmp20 = data((data(:,1)>0.0390),:);
tmp20 = tmp20((tmp20(:,1)<0.0392),:);
tmp20 = tmp20((tmp20(:,2)==2),:);
tmp30 = data((data(:,1)>0.0507),:);
tmp30 = tmp30((tmp30(:,2)==2),:);

figure
hold on
plotthis=tmp05;
plot(plotthis(:,3),plotthis(:,1)-plotthis(:,1),'k--');
plot(plotthis(:,3), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'+');
plot(plotthis(:,3),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'o')

plotthis=tmp20;
plot(plotthis(:,3), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'x');
plot(plotthis(:,3),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'s')

plotthis=tmp30;
plot(plotthis(:,3), (plotthis(:,4)-plotthis(:,1))./plotthis(:,1),'d');
plot(plotthis(:,3),(polyval(p_gi,plotthis(:,6))-plotthis(:,1))./plotthis(:,1),'v')

legend('ideal','parametric 50 ml','GI 50 ml','parametric 250 ml','GI 250 ml','parametric 550 ml','GI 550 ml');
title('parametric vs. GI');
ylabel('relative error')
xlabel('SNR');