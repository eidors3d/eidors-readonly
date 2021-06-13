function s_mat = system_mat_instrument( img);
% SYSTEM_MAT_INSTRUMENT: 
%  Calculate systems matrix inclusing modelling of instrument impedance
% img.fwd_model.system_mat = @system_mat_instrument
% img.fwd_model.system_mat_instrument.clist = [CONNECTION LIST]
% 
% where
%  connection list =
%     [node1a, node1b, admittance_1ab]
%     [node2a, node2b, admittance_2ab]

% (C) 2021 Andy Adler. License: GPL version 2 or version 3
% $Id$
   new_c_list = img.fwd_model.system_mat_instrument.connect_list;

   s_mat= system_mat_1st_order(img);
   E = s_mat.E;

   n_nodes = num_nodes(img);
   if ~isempty(new_c_list);
      n_elecs = max([max(new_c_list(:,1:2)),num_elecs(img)]);
   else
      n_elecs = num_elecs(img);
   end
   n_max  = n_nodes + n_elecs; 
   Eo = sparse(n_max, n_max);
   Eo(1:size(E,1),1:size(E,2)) = E;
   for i=1:size(new_c_list,1);
     x = new_c_list(i,1) + n_nodes;
     y = new_c_list(i,2) + n_nodes;
     c = new_c_list(i,3);
     Eo(x,y) = Eo(x,y) - c;
     Eo(y,x) = Eo(y,x) - c;
     Eo(x,x) = Eo(x,x) + c;
     Eo(y,y) = Eo(y,y) + c;
   end
   
   s_mat.E = Eo;
