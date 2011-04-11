function [snr,accuracy,sym_err]=data_analysis(Reshape)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [snr,accuracy,sym_err]=data_analysis(Reshape)
% Name:data_analysis.m
% input: Reshape (0/1) is a status parameter used to check if reshape of 
% signal from sin to plateau data shape is needed (1) or not (0)
%
% output: snr, accuracy,sym_err
% sym_err: calculate symmetrical measurement electrode errors

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

% example: [snr, accuracy,sym_err]= data_analysis(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load measurement and model data
Vmeas = load_meas(Reshape);
[Vsim, fmdl] = model_data;

%% calculate snr
snr = 20*log10(abs(mean(Vmeas,2)./ std(Vmeas,0,2)));
snr_max=max(snr); snr_min=min(snr); snr_mean=mean(snr);

%%  calculate accuracy
accuracy = calc_accuracy(Vmeas, Vsim);

%% calculate symmetrical measurement error
[sym_err, tot_sym_err] = calc_sym_err(Vmeas);

%% Plotting
plott(snr,accuracy,sym_err);
end

function accuracy = calc_accuracy(Vmeas, Vsim)

Vm = mean(Vmeas,2);
% normalized to 0 to 1
Vm = (Vm(:) - min(Vm)) /(max(Vm)-min(Vm));
% normalized to zero to 1
Vs = (Vsim.meas(:) - min(Vsim.meas))/(max(Vsim.meas)-min(Vsim.meas));
% Accuracy
accuracy = (1-abs(Vm - Vs))*100;
end

function Vmeas = load_meas(Reshape)
if Reshape
    % if data is in the Sheffield measurement sequence
    datapath =  ('.././data/path/');
    vh = load([datapath,'homo.mat']);
    Vh =abs([vh.Eit_Data{:}]);
    Vmeas= sin2plateau(V)
else  % if data is plateau data shape 
    datapath =  ('.././data/accuracy/');
    vh = load([datapath,'post0.mat']);
    Vh =abs([vh.Eit_Data{:}]);
    IDX = [18:2:416,2:2:16];
    Vmeas= Vh(IDX,:);
end

end

%%
function [Vsim, fmdl] = model_data
% simulate measurement voltageabs
h= 30/2; w = 28;
curr = 1; % applied current in mA
fmdl = ng_mk_cyl_models([2*h,w/2, 0.6],[16,h],[0.05]);
%save Tank_model.mat  fmdl
% load('.././data/accuracy/Tank_model.mat')

%  Stimulation patters
k=0; 
for i=2:14
    for j=1:16; k=k+1;
        fmdl.stimulation(k).stim_pattern=sparse(rem(j+[-1;0],16)+1,1,[1;-1],16,1);
        fmdl.stimulation(k).meas_pattern=sparse(rem(i+j+[-1;0],16)+1,1,[-1;1],16,1)';
    end
end
img = mk_image(fmdl,1);
Vsim = fwd_solve(img);

end
%%
function [sym_err, tot_sym_err] = calc_sym_err(v)
%Procedures for the code 
%1.	Average multiple frames for a complete set of 208 measurements
%2.	Sort symmetric measurements into two groups
%3.	Take the RMS from two group of symmetric voltages

v = mean(v,2);
nmeas=size(v,1);

%rearrange data sequence
count=0;
for j=1:13
    for i=1:16
        count =1+count;
        vv(13*(i-1)+j) = v(count);
    end
end

for i=1:nmeas
    if mod(i,13) ~= 0
        idx(i) = mod(i,13);
    else
        idx(i) = 13;
    end
end

count1 = 0; count2 = 0;
for j = 1:nmeas
    if mod(idx(j),7) ~=0
        
        if idx(j) <= 6
            count1 = count1+1;
            VR(count1) = vv(j);
        else
            count2 = count2+1;
            VL1(count2) = vv(j);
        end
        
    end
end

temp(1:6)=0;
for k = 1:6:count2
    temp(1:6) = VL1(k:k+5);
    VL(k:k+5) = temp(6:-1:1);
end

% plot(VR,'DisplayName','VR','YDataSource','VR');
% hold all;
% plot(VL,'DisplayName','VL','YDataSource','VL');
% hold all;;hold off;figure(gcf);
sym_sum = 0;
for i = 1:count2
    sym_err(i) = ((VL(i)- VR(i))/VR(i));
    sym_sum = sym_sum + (sym_err(i))^2;
end
tot_sym_err = sqrt(sym_sum/count2)
end


function [symer, tot_symer] = calc_sym_err2(v)
% Approach 2: not fully tested
% symmetricity with respect to the middle of the x-axis
% the signal is in plateau signal shape
symm = 0;
for ii = 1:104
    symer (ii) = (v(ii) - v(208-ii))/v(ii);
    symm = symm + (symer(ii))^2;
end
tot_symer = sqrt(symm/104);
end

%%
function plott(snr,accuracy,sym_err)
figure;
subplot(1,3,1);
plot(snr,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',2);
ylabel('SNR values (dB)');
xlabel('Channel number');
axis( [ 0 215 40 80 ] );

subplot(1,3,2);
plot(accuracy,'-bs','LineWidth',1,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',1.5);
ylabel('Accuracy');
xlabel('Channel number');
axis( [ 0 215 88 100 ] );

%% plotting symmetric electrode measurement error
subplot(1,3,3);
plot(sym_err,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w',...
    'MarkerSize',4);
xlabel('Symmetric measurements');
ylabel(' Symmetric measurement error');

end
