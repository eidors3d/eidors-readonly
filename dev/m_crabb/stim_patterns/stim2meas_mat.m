function [meas_mat]=stim2meas_mat(stim)
%Function taking eidors stimulation pattern structure to measurement matrix
%stim - Eidors stimulation structure
%meas_mat - Each row is indices corresponding to [+I,-I,+V,-V]
%nelecs - total number of electrodes
%NOTE - The function only allows a +I,-I stimulation without patches

nstims=size(stim,2); idx=0;
for i=1:nstims   
   %Get the stimulation nodes
   stim_pat=stim(i).stim_pattern;
   %Get the measurement patterns and number of measurements
   meas_pat=stim(i).meas_pattern;  n_meas_i=size(meas_pat,1);
   
   %Find the positive/negative current indices and assign
   posI=find(stim_pat(:) == 1); meas_mat(idx+(1:n_meas_i),1)=posI;
   negI=find(stim_pat(:) ==-1); meas_mat(idx+(1:n_meas_i),2)=negI;
     
   %Loop through measurements and assign 3rd/4th column
   for j=1:n_meas_i
      posV=find(meas_pat(j,:) == 1); meas_mat(idx+j,3)=posV;
      negV=find(meas_pat(j,:) ==-1); meas_mat(idx+j,4)=negV;
   end
   idx=idx+n_meas_i;    
end

end