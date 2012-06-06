function img= nodal_solve( inv_model, data1, data2)
% AA_INV_SOLVE inverse solver using approach of Adler&Guardo 1996
% img= nodal_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$


dv = calc_difference_data( data1, data2, inv_model.fwd_model);

sol = get_RM( inv_model ) * dv;

% create a data structure to return
img.name= 'solved by nodal_solve';
img.node_data = sol;
img.fwd_model= inv_model.fwd_model;

function one_step_inv = get_RM( inv_model );
   fwd_model= inv_model.fwd_model;

   % The one_step reconstruction matrix is cached
   one_step_inv = eidors_obj('get-cache', inv_model, 'nodal_solve');
   if ~isempty(one_step_inv)
       eidors_msg('nodal_solve: using cached value', 3);
       return;
   end

   img_bkgnd= calc_jacobian_bkgnd( inv_model );
   J = calc_jacobian( fwd_model, img_bkgnd);

   RtR = calc_RtR_prior( inv_model );
   W   = calc_meas_icov( inv_model );
   hp  = calc_hyperparameter( inv_model );

   [e2n,Ne, Nn] = elem2node( fwd_model.elems );
   if size(J,2) == Ne;
      J= J*e2n;
   end
   if size(RtR,2) == Ne;
      RtR = e2n'*RtR*e2n;
      RtR = RtR +  1e-4*spdiags(diag(RtR),0,Nn,Nn); %Need just a little regularization
   end

   one_step_inv= (J'*W*J +  hp^2*RtR)\J'*W;

   eidors_obj('set-cache', inv_model, 'nodal_solve', one_step_inv);
   eidors_msg('nodal_solve: setting cached value', 3);

function [e2n, Ne, Nn] = elem2node( elems )
   [Ne,d] = size(elems);
   elemo = (1:Ne)'*ones(1,d);
   e2n = sparse(elemo,elems,1/d); 
   Nn = size(e2n,2);
