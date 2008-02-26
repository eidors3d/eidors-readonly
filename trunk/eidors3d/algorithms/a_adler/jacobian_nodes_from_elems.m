function J= jacobian_nodes_from_elems( fwd_model, img);
% JACOBIAN_NODES_FROM_ELEMS: calculate Jacobian on Nodes from Elem solver
% Calculate Jacobian Matrix for EIT Alg of Adler & Guardo 1996
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: jacobian_nodes_from_elems.m,v 1.1 2008-02-26 21:23:16 aadler Exp $

EtoN = mapper_elems_nodes( fwd_model);
n_nodes= size(EtoN,1);

iEtoN= (EtoN'/(EtoN*EtoN'+1e-6*speye(n_nodes)));

% Create an image on the elements with a fwd_model on the elemtents
img_e= img;
img_e.elem_data= 


return;
% Three different ways to invert iEtoN. We only need something fairly
% approximate.
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
