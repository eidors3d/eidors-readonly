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



save(savefile,'test_v');


% Radial Movement - $Id$  
function [xp,yp,zp]= simulation_radmove(f_frac, radius, z0,zt);
   rp= f_frac*radius; 
%  cv= 2*pi*f_frac * 73;
   cv= 0.5;
   xp= rp * cos(cv);
   yp= rp * sin(cv);
   if nargin==4; zp = mean([zt,z0]); else; zp= 0; end


