function J= jacobian_nodes_from_elems( fwd_model, img);
% JACOBIAN_NODES_FROM_ELEMS: calculate Jacobian on Nodes from Elem solver
% Calculate Jacobian Matrix for EIT Alg of Adler & Guardo 1996
% J         = Jacobian matrix
% fwd_model = forward model defined on nodes (elems may not be defined)
% jacobian_nodes_from_elems.fwd_model 
%           = full forward model (with nodes and elements)
%
% img = image background for jacobian calc

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: jacobian_nodes_from_elems.m,v 1.2 2008-02-26 22:05:02 aadler Exp $

fem_fmdl= fwd_model.jacobian_nodes_from_elems.fwd_model;
EtoN = mapper_nodes_elems( fem_fmdl);
n_elems= size(EtoN,2);

% Create an image on the elements with a fwd_model on the elemtents
img_e= img;
% jnkflag needed to make lsqr shut up
mu=1e-3;
[e_d, jnkflag] = lsqr([EtoN;mu*speye(n_elems)], ...
                      [img_e.node_data;zeros(n_elems,1)], ...
                       1e-8,1000); % should be enough, accuracy is ok
img_e.elem_data= e_d;
img_e.fwd_model = fem_fmdl;

J= calc_jacobian(fem_fmdl, img_e)*EtoN';

return;
% Three different ways to invert iEtoN. We only need something fairly
% approximate. This is test code to check
[Nn,Ne]= size(EtoN);  
n1= ones(Nn,1); mu=1e-3;
t=cputime;
   ed1= (EtoN'/(EtoN*EtoN'+mu^2*speye(Nn)))*n1;
disp([cputime-t, std(ed1)]);
t=cputime;
   ed2= ((EtoN'*EtoN+mu^2*speye(Ne))\EtoN')*n1;
disp([cputime-t, std(ed2)]);
t=cputime;
% Use model saying [E;m*I]*i= [n;0]   
   ed3= lsqr([EtoN;mu*speye(Ne)],[n1;zeros(Ne,1)],1e-8,1000);
disp([cputime-t, std(ed3)]);
