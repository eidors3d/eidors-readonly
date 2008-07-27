% Human breathing $Id$

% Electrodes on back
load montreal_data_1995
v(4).vh = double( zc_resp(idx,1) );
v(4).vi = double( zc_resp(idx,22) );
