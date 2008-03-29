% Create animated graphics $Id: simulate_move2_05.m,v 1.1 2008-03-29 15:39:56 aadler Exp $

% Trim images
!find -name 'simulate_move2_04a*.png' -exec convert  -trim '{}' PNG8:'{}' ';'

% Convert to animated Gif
!convert -delay 50 simulate_move2_04a*.png -loop 0 simulate_move2_05a.gif
