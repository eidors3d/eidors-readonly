% create model $Id: simulate_move2_01.m,v 1.2 2008-03-29 14:59:19 aadler Exp $

% Model parameters
n_elec= 16;
n_nodes= 200;

% Create electrodes
refine_level=4; %electrode refinement level
elec_width= .1;
z_contact = 0.01;
% elect positions
th=linspace(0,2*pi,n_elec(1)+1)';th(end)=[];
elec_posn= [sin(th),cos(th)];
[elec_nodes, refine_nodes] = dm_mk_elec_nodes( elec_posn, ...
       elec_width, refine_level);

% Define circular medium
fd=inline('sum(p.^2,2)-1','p');
bbox = [-1,-1;1,1];
smdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, ...
                          elec_nodes, refine_nodes, z_contact);
