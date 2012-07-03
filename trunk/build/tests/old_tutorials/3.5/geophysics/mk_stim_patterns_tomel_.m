function [st,data]= mk_stim_patterns_tomel(elec_posn,data_tomel,name,configuration)
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
% Dominique Gibert, June 2007
%
%%
A= data_tomel(:,3); % Positive current electrode number
B= data_tomel(:,4); % Negative current electrode number
M= data_tomel(:,5); % Positive voltage electrode number
N= data_tomel(:,6); % Negative voltage electrode number
current= data_tomel(:,7)/1000.; % injected current (mA -> A)
E= elec_posn(:,1);  % List of electrodes in the same order as in the
%construction of the cross-section of the gallery
AB_all= A+i*B;      % Use complex number to simultaneously process A and B electrodes
MN_all= M+i*N;
AB_uni= unique(AB_all); % Keep only one copy of each (A,B) configuration
AB_uni= AB_uni(ismember(real(AB_uni),E)); % Keep only those A present in the E list
AB_uni= AB_uni(ismember(imag(AB_uni),E)); % Keep only those B present in the E list
n_stim= length(AB_uni); % OK, the AB list is now minimal and coherent with the E list
A_uni= real(AB_uni);
B_uni= imag(AB_uni);
n_elec= length(E);
place_data= [];
for k = 1:n_stim
    stim_pat= sparse(n_elec,1);
    [tf,loc]= ismember(A_uni(k),E);
    stim_pat(loc)= +current(k); % current is not the right one !!!
    [tf,loc]= ismember(B_uni(k),E);
    stim_pat(loc)= -current(k);
    meas_for_ABk= find(AB_all == AB_uni(k));
    Mk= M(AB_all == AB_uni(k));
    Nk= N(AB_all == AB_uni(k));
    placeABk= find(AB_uni(k)==AB_all);
    n_meas= size(meas_for_ABk,1);
    meas_pat= sparse(n_meas,n_elec);
    for l = 1:n_meas
        placeMNk= find(Mk(l)+i*Nk(l)==MN_all);
        place_data= [place_data; intersect(placeMNk,placeABk)];
        [tf,loc]= ismember(Mk(l),E);
        meas_pat(l,loc)= +1;
        [tf,loc]= ismember(Nk(l),E);
        meas_pat(l,loc)= -1;
    end
    st(k).stimulation= 'mA';
    st(k).delta_time= 0.0;
    st(k).stim_pattern= stim_pat;
    st(k).meas_pattern= meas_pat;
end
% Now re-order the data list according to the stimulation pattern
data.meas= data_tomel(place_data,8);
data.name= name;
data.type= 'data';
data.time= -1;
data.configuration= configuration;
end

% function st= mk_stim_patterns_tomel(elec_posn,data_tomel)
% %% Returns the stimulation sub-structure corresponding to a tomel
% %  measurement file.
% %
% % elec_posn = 4-column matrix as for example:
% %
% %    EZG04_Ring1= [ ...
% %         101   1.037202       0.000871      -0.465710  ... 
% %         132   1.000000      -0.650000      -0.440000];    
% %
% % where:
% %   col 1 = electrode number (arbitrary)
% %   col 2 = x position of the electrode (m) (not used)
% %   col 3 = y position of the electrode (m) (not used)
% %   col 4 = z position of the electrode (m) (not used)
% %
% % Only columns 1 (i.e. the electrode number) is used in this function
% % If called outside of mk_gallery, elec_posn must be the identical to the
% % one sent in mk_gallery when building the gallery FEM structure
% %
% % data_tomel = TOMEL data are given as a matrix with at least 7 colums as,
% % for example:
% %
% %  Data_Ring1_Dec2004_Wen32_1= [ ...
% %       1     223    101    104    102    103   100    0.260163 ...
% %     160    2285    132    115    105    110   100    0.041110];
% %
% % where:
% %   col 1 = measurement number (not used)
% %   col 2 = time elapsed since begining of measurement sequence (in seconds, not used)
% %   col 3 = A electrode ( = positive current electrode)
% %   col 4 = B electrode ( = negative current electrode)
% %   col 5 = M electrode ( = positive potential electrode)
% %   col 6 = N electrode ( = negative potential electrode)
% %   col 7 = intensity (mA) of the current injected from A to B
% %   col 8 = voltage measured V(M) - V(N) for the measurement pattern (not used)
% %
% % In this function, the intensities are converted from mA to A
% %
% % Dominique Gibert, April 2007
% %
% %%
% A= data_tomel(:,3); % Positive current electrode number
% B= data_tomel(:,4); % Negative current electrode number
% M= data_tomel(:,5); % Positive voltage electrode number
% N= data_tomel(:,6); % Negative voltage electrode number
% current= data_tomel(:,7)/1000.; % injected current (mA -> A)
% E= elec_posn(:,1); % List of electrodes in the same order as in the
% %construction of the cross-section of the gallery
% n_stim= length(A);
% n_elec= length(E);
% for i = 1:n_stim
%     stim_pat= zeros(n_elec,1);
%     [tf,loc]= ismember(A(i),elec_posn(:,1));
%     stim_pat(loc)= +current(i);
%     [tf,loc]= ismember(B(i),elec_posn(:,1));
%     stim_pat(loc)= -current(i);
%     meas_pat= zeros(n_elec,1);
%     [tf,loc]= ismember(M(i),elec_posn(:,1));
%     meas_pat(loc)= +1;
%     [tf,loc]= ismember(N(i),elec_posn(:,1));
%     meas_pat(loc)= -1;
%     st(i).stimulation= 'mA';
%     st(i).delta_time= 0.0;
%     st(i).stim_pattern= stim_pat;
%     st(i).meas_pattern= meas_pat';
% end
% end
