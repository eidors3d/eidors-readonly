function st= mk_stim_patterns_tomel(elec_posn,data_tomel)
%% Returns the stimulation sub-structure corresponding to a tomel
%  measurement file.
%
% elec_posn = 4-column matrix as for example:
%
%    EZG04_Ring1= [ ...
%         101   1.037202       0.000871      -0.465710  ... 
%         132   1.000000      -0.650000      -0.440000];    
%
% where:
%   col 1 = electrode number (arbitrary)
%   col 2 = x position of the electrode (m) (not used)
%   col 3 = y position of the electrode (m) (not used)
%   col 4 = z position of the electrode (m) (not used)
%
% Only columns 1 (i.e. the electrode number) is used in this function
% If called outside of mk_gallery, elec_posn must be the identical to the
% one sent in mk_gallery when building the gallery FEM structure
%
% data_tomel = TOMEL data are given as a matrix with at least 7 colums as,
% for example:
%
%  Data_Ring1_Dec2004_Wen32_1= [ ...
%       1     223    101    104    102    103   100    0.260163 ...
%     160    2285    132    115    105    110   100    0.041110];
%
% where:
%   col 1 = measurement number (not used)
%   col 2 = time elapsed since begining of measurement sequence (in seconds, not used)
%   col 3 = A electrode ( = positive current electrode)
%   col 4 = B electrode ( = negative current electrode)
%   col 5 = M electrode ( = positive potential electrode)
%   col 6 = N electrode ( = negative potential electrode)
%   col 7 = intensity (mA) of the current injected from A to B
%   col 8 = voltage measured V(M) - V(N) for the measurement pattern (not used)
%
% In this function, the intensities are converted from mA to A
%
% Dominique Gibert, April 2007
%
%%
A= data_tomel(:,3); % Positive current electrode number
B= data_tomel(:,4); % Negative current electrode number
M= data_tomel(:,5); % Positive voltage electrode number
N= data_tomel(:,6); % Negative voltage electrode number
current= data_tomel(:,7)/1000.; % injected current (mA -> A)
E= elec_posn(:,1); % List of electrodes in the same order as in the
%construction of the cross-section of the gallery
n_stim= length(A);
n_elec= length(E);
for i = 1:n_stim
    stim_pat= zeros(n_elec,1);
    [tf,loc]= ismember(A(i),elec_posn(:,1));
    stim_pat(loc)= +current(i);
    [tf,loc]= ismember(B(i),elec_posn(:,1));
    stim_pat(loc)= -current(i);
    meas_pat= zeros(n_elec,1);
    [tf,loc]= ismember(M(i),elec_posn(:,1));
    meas_pat(loc)= +1;
    [tf,loc]= ismember(N(i),elec_posn(:,1));
    meas_pat(loc)= -1;
    st(i).stimulation= 'mA';
    st(i).delta_time= 0.0;
    st(i).stim_pattern= stim_pat;
    st(i).meas_pattern= meas_pat';
end
return;