Original = eidors_readdata('nm18.get','GET-RAW');
Orig_ref = mean(Original,2); % reference frame - average of all frames
Orig_imag=inv_solve(inv_mdl,Original(:,691:717),Orig_ref); %reconstruct difference images

%Corrected = eidors_readdata('nm18_corrected.get','GET-RAW');  % If using GUI
Corrected = frame_timing_correct( Original, 'FFT');
Corr_ref = mean(Corrected,2); % reference frame - average of all frames
Corr_imag=inv_solve(inv_mdl,Corrected(:,691:717),Corr_ref); %reconstruct difference images
