% Create and show square model $Id: dm_geophys02.m,v 1.1 2008-04-03 19:33:28 aadler Exp $

% Create square mesh model
imdl= mk_common_model('c2s',16);
cmdl= rmfield(imdl.fwd_model,{'electrode','stimulation'});

% magnify and move down onto geophysics model
cmdl.nodes(:,2)= cmdl.nodes(:,2) - 1.05;
cmdl.nodes= cmdl.nodes*4;

% assign one parameter to each square
e= size(cmdl.elems,1);
params= ceil(( 1:e )/2);
cmdl.coarse2fine = sparse(1:e,params,1,e,max(params));

clf;
show_fem(fmdl);
hold on;
h= trimesh(cmdl.elems,cmdl.nodes(:,1),cmdl.nodes(:,2));
set(h,'Color',[0,0,1],'LineWidth',2);
hold off
axis image

print -dpng -r125 dm_geophys02a.png
