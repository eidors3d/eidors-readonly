% create model $Id: simulate_move2_01.m,v 1.1 2008-03-29 03:03:14 aadler Exp $

n_elec= 16;
n_nodes= 200;
refine_level=2; %electrode refinement
elec_width= .1;
th=linspace(0,2*pi,n_elec(1)+1)';th(end)=[];
elec_posn= [sin(th),cos(th)];
[elec_nodes, refine_nodes] = dm_mk_elec_nodes( elec_posn, ...
       elec_width, refine_level);

th=linspace(0,2*pi,20)'; th(end)=[];
%fd=inline('max(sum(p.^2,2)-1,0.098-sqrt(sum(p.^2,2)))','p');
fd=inline('sum(p.^2,2)-1','p');
bbox = [-1,-1;1,1];
z_contact = 0.01;
stim_pat= mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},1);
smdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, ...
                          elec_nodes, refine_nodes, z_contact);

smdl.stimulation= stim_pat;
img1= eidors_obj('image','','fwd_model',smdl,...
                 'elem_data',ones(size(smdl.elems,1),1) );

refine_nodes= [refine_nodes;0.10*[sin(th)+1,cos(th)-1]];
smdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, ...
                          elec_nodes, refine_nodes, z_contact);

smdl.stimulation= stim_pat;
img2= eidors_obj('image','','fwd_model',smdl,...
                 'elem_data',ones(size(smdl.elems,1),1) );

show_fem(smdl)

