function img= aa_inv_solve( inv_model, data1, data2)
% BACKPROJ_SOLVE inverse solver using backprojection
% img= backproj_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% inv_model.backproj_solve.type = 'naive' (DEFAULT)
%    use naive (unfiltered algorithm)
% inv_model.backproj_solve.type = 'filtered' (NOT IMPLEMENTED YET)
%    ref: Barber DC Brown BH, "fast reconstruction of resistance
%         images", clin Phys Physiol Mes, pp 47-54, vol 8,sup A,1987
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: backproj_solve.m,v 1.1 2007-10-22 20:04:07 aadler Exp $

try
   type= inv_model.backproj_solve.type;
catch
   type= 'naive';
end

fwd_model= inv_model.fwd_model;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'backproj_solve');
if ~isempty(one_step_inv)
    eidors_msg('backproj_solve: using cached value', 2);
else
    Jbp = calc_backprojection_mask( fwd_model );

    one_step_inv= (J'*W*J +  hp^2*RtR)\J'*W;

    eidors_obj('set-cache', inv_model, 'backproj_solve', one_step_inv);
    eidors_msg('backproj_solve: setting cached value', 2);
end

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

sol = one_step_inv * dv;

% create a data structure to return
img.name= 'solved by backproj_solve';
img.elem_data = sol;
img.fwd_model= fwd_model;

function Jbp = calc_backprojection_mask( fmdl );
   % create homog image
   himg= eidors_obj('image','', 'fwd_model',fmdl, ...
                    'elem_data',ones(size(fmdl.elems,1),1) );
   node_v= calc_all_node_voltages( himg );

   pp= aa_fwd_parameters( fmdl );

   elem_v= reshape(node_v( fmdl.elems',:),pp.n_dims+1,[],pp.n_elec);
   elem_v= squeeze(mean(elem_v,1));

   meas_v= pp.N2E * node_v;

   stimulation= fwd_model.stimulation;
   for i = 1:length(stimulation);
     [idx_pl,jnk]=find(fmdl.stimulation(1).meas_pattern'<0)
     [idx_mi,jnk]=find(fmdl.stimulation(1).meas_pattern'>0)


keyboard
   sel= 4;
   ed= elem_v(:,sel);
   img.elem_data= ed;
   subplot(2,3,idx);
   show_fem(img);

   ej = zeros(size(ed));
   for i=1:16
     ej= ej + ( ed < meas_v(i,sel) );
   end
   img.elem_data= ej;
   subplot(2,3,idx+3);
   show_fem(img);
e

function N2E= calc_mappings(fwd_model);
   electrode= fwd_model.electrode;
   n_elec= length(electrode);
   N2E = sparse(n_elec, size(fwd_model.nodes,1));
   for i=1:n_elec
       elec_nodes = electrode(i).nodes;
       if length(elec_nodes) == 1 
          N2E(i, elec_nodes) = 1;
       else
          srf_area   = get_srf_area( bdy, elec_nodes, fwd_model.nodes);
          N2E(i, elec_nodes) = srf_area/length(elec_nodes);
       end
   end
