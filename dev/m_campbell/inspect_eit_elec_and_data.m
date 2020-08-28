function inspect_eit_elec_and_data(eit_file, imdl, thresh)
% -------------------------------------------------------------------------
% Description:
%   inspect_eit_elec_and_data(eit_file, imdl, thresh)
% -------------------------------------------------------------------------
% Parameters:
%   eit_file:
%       Target data file ending in .eit file extension or a cell containing
%       the outpus of eidors_readdata {data, auxdata}.
%   imdl:
%       EIDORS inverse model structure
%   thresh:
%       Electrode score at which to label an electrode poor.
% ------------------------------------------------------------------------- 
% Returns:
%   Total boundary voltage plot
%   IQ plot
%   U-shape plot
%   Frequency spectrum plot
%   Contact impedance plot
%   Average measurement for measurement pair plot
%   Highlights low-quality electrodes and measurements
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
msel= imdl.fwd_model.meas_select;
mm  = find(msel);

if nargin == 2
   thresh = 0.25; 
end % end if

if isstruct(eit_file) && strcmp(thresh, 'hamburg')
    data            = eit_file.data;
    elec_impedance  = eit_file.elec_impedance;
    fs              = eit_file.fs;
    thresh          = 0.25;
elseif iscell(eit_file) % input is output of eidors_read_data
    if length(eit_file)==2
        if isstruct(eit_file{2})
            data    = eit_file{1};
            auxdata = eit_file{2};
        else
            auxdata = eit_file{1};
            data    = eit_file{2};
        end % end if
        fs              = 1e6./median(diff(auxdata.t_rel));     % framerate is median dif of time points/ 1000000 (convert to s)
        impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
        elec_impedance  = auxdata.elec_impedance* impedanceFactor;
    else
        disp("Whoops!");
    end % end if
elseif strcmp(eit_file(end-3:end), '.eit')
    [data, auxdata] = eidors_readdata(eit_file);
    fs              = 1e6./median(diff(auxdata.t_rel));         % framerate is median dif of time points/ 1000000 (convert to s)
    impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15;   % = 0.9633 / 2^15;
    elec_impedance  = auxdata.elec_impedance* impedanceFactor;
else
    data = eit_file;
end % end if

n_samples       = size(data, 2);
data_samp       = round(1:(n_samples/100):n_samples);
elec_impedance  = abs(elec_impedance);

[wElecs, scores, mmScores, ~] = worst_n_elecs(data, imdl);
% add plot bad measurements
badMeas     = mmScores > thresh;
badElecIdx  = scores > thresh;
if length(badElecIdx) >= 1
    badElecs= wElecs(badElecIdx);
    kk      = meas_icov_rm_elecs(imdl, badElecs);
    ee1     = find(diag(kk)~=1);
end % end if
ee = badMeas;
ei = mean(elec_impedance, 2);

% total boundary voltage
subplot(4,2,1:2);
    vv = detrend(real(data));
    xax= (1:size(vv, 2))/ fs;
    plot(xax, sum(vv(msel,:), 1));
    ylabel('Voltage (V)')
    xlabel('Time (s)');
    title('Total Boundary Voltage (Raw)');

% U shapes 
subplot(424);
    vk = 1e3 * mean(vv, 2);
    plot(vk,'k');    % all meas
    hold on;
    vk(~msel,:) = NaN;
    plot(vk,'b');   % selected meas
    vn = NaN * vk; vn(ee1,:) = vk(ee1,:); % measurements from rejected electrodes
    plot(vn,'mo');
    vn = NaN * vk; vn(ee,:) = vk(ee,:); % independently rejected measurements
    plot(vn,'ro');
    hold off;
    box off; 
    xlim([0,400]);
    title 'U shapes';

% IQ plot
subplot(4,2,[3,5,7]);
    % all meas
    plot(1e3 * data(:, data_samp), 'k+');
    hold on;
    % used meas
    plot(1e3 * data(msel, data_samp), 'b+');
        
    ee = mm(ee);
    plot(1e3 * data(ee,data_samp), 'r+');   % independently rejected measurements
    
    ee1= mm(ee1);
    plot(1e3 * data(ee1,data_samp), 'm+');  % measurements from rejected electrodes (covers rejected meas)

    hold off; 
    box off;
    axis equal;
    title 'IQ plot';

% median contact impedance
subplot(426);
    b           = bar(ei,'FaceColor','flat');
    b.CData(:,:)= repmat([0 0 1], 32, 1);
    be          = badElecs;
    for e= 1:length(be)
        b.CData(be(e), :)= [1 0 0];
    end
    hold on;
    eb = errorbar(1:32, ei, max(elec_impedance, [], 2)-ei, min(elec_impedance, [], 2)- ei);
    set(eb,'Color',[0, 0, 0],'LineStyle','none');
    hold off;
    box off; 
    xlim([0, 33]);
    title 'Median Elec Z';

% Fourier Series
subplot(428);
    ll  = size(data, 2);    
    ft  = fft(detrend(data.').', [], 2);
    fax = linspace(0, fs, ll+1); fax(end)=[];
    semilogy(fax, mean(abs(ft)));
    ylim([1e-4, 1]);
    box off; 
    xlim([0, fs/2]);
    title 'Fourier Spectrum (Hz)';

end % end function

