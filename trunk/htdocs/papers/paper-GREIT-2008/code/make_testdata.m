function make_testdata( savefile )
% Make test data and save in file
% 
% (C) Andy Adler 2008 $Id$
% Licensed under GNU GPL v2 or v3

% exclude measures at electrodes
[x,y]= meshgrid(1:16,1:16); idx= abs(x-y)>1 & abs(x-y)<15;

K=1;

% PHANTOM DATA
load iirc_data_2006
test_v(K).vh= - real(v_reference(idx,1));
test_v(K).vi= - real(v_rotate(idx,1));

load montreal_data_1995
K=K+1;
test_v(K).vh= double( zc_h_demo4(idx) );
test_v(K).vi= double( zc_demo4(idx,1) );
K=K+1;
test_v(K).vh= double( zc_h_demo4(idx) );
test_v(K).vi= double( zc_demo4(idx,5) );
K=K+1;
test_v(K).vh= double( zc_h_demo4(idx) );
test_v(K).vi= double( zc_demo4(idx,10) );


%  HUMAN TIDAL BREATHING

% Electrodes on back

K=K+1;
test_v(K).vh = double( zc_resp(idx,1) );
test_v(K).vi = double( zc_resp(idx,22) );

% Contrib EIT data - PIG
vv= eidors_readdata('p1130107.get');
vh = mean(vv(idx,1:40),2);
K=K+1;
test_v(K).vh = vh;
%test_v(K).vi = vv(idx,1345); % PEEP 5 end-expi
test_v(K).vi = vv(idx,387); % PEEP 10 end-expi incr
K=K+1;
test_v(K).vh = vh;
%test_v(K).vi = vv(idx,1145); % PEEP 15 end-expi
test_v(K).vi = vv(idx,1257); % PEEP 10 end-expi decr

% Contrib EIT data - PIG
vv= eidors_readdata('1-control.RAW');
vh= mean(vv(idx,[37, 42, 46, 50, 54, 58, 62, 66, 70, 75]),2);
vi= mean(vv(idx,[35, 39, 43, 47, 51, 64, 68, 72, 76, 80]),2);
K=K+1;
test_v(K).vh = vh; % FRC
test_v(K).vi = vi; % VT

vv= eidors_readdata('1-injury.RAW');
vh= mean(vv(idx,[37, 41, 45, 58, 62, 66, 70, 74, 82, 87]),2);
vi= mean(vv(idx,[34, 47, 51, 55, 59, 63, 76, 80, 84, 88]),2);
K=K+1;
test_v(K).vh = vh; % FRC
test_v(K).vi = vi; % VT

% Contrib EIT data - Pneumothroax / Pleural Effusion
vv= eidors_readdata('goev354005.get');
vh= mean(vv(idx,[21,60,97,134,206,241,278,316,392,430]),2);

vv= eidors_readdata('goev354003_512-1560.get');
va= mean(vv(idx,[569,605,643,680,715,753,827,866,901,940]),2);

vv= eidors_readdata('goev354008_947-1559.get');
vf= mean(vv(idx,[51, 89, 125, 163, 201, 236, 274, 311, 349, 385, 423, 461, 496, 534]),2);

K=K+1;
test_v(K).vh = vh; % FRC
test_v(K).vi = va; % Air - FRC
K=K+1;
test_v(K).vh = vh; % FRC
test_v(K).vi = vf; % Fluid - FRc



save(savefile,'test_v');

