function img= inv_solve_backproj( inv_model, data1, data2)
% INV_SOLVE_BACKPROJ inverse solver using backprojection
% NOTE: This is the beginnings of an attempt to reproduce
%  the backprojection algorithm implemented in the
%  Sheffield MKI EIT system. It is far from complete.
%
% If you wish to use the actual algorithm, use the
%  function "mk_common_gridmdl('backproj')"
%
% img= inv_solve_backproj( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% inv_model.inv_solve_backproj.type = 'naive' (DEFAULT)
%    use naive (unfiltered algorithm)
% inv_model.inv_solve_backproj.type = 'filtered' (NOT IMPLEMENTED YET)
%    ref: Barber DC Brown BH, "fast reconstruction of resistance
%         images", clin Phys Physiol Mes, pp 47-54, vol 8,sup A,1987
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id$

try
   type= inv_model.inv_solve_backproj.type;
catch
   type= 'naive';
end

fwd_model= inv_model.fwd_model;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'inv_solve_backproj');
if ~isempty(one_step_inv)
    eidors_msg('inv_solve_backproj: using cached value', 2);
else
    Jbp = calc_backprojection_mask( fwd_model, type );

    one_step_inv= Jbp;

    eidors_obj('set-cache', inv_model, 'inv_solve_backproj', one_step_inv);
    eidors_msg('inv_solve_backproj: setting cached value', 2);
end

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

sol = one_step_inv * dv;

% create a data structure to return
img.name= 'solved by inv_solve_backproj';
img.elem_data = sol;
img.fwd_model= fwd_model;

function Jbp = calc_backprojection_mask( fmdl , type);
   % create homog image
   himg= eidors_obj('image','', 'fwd_model',fmdl, ...
                    'elem_data',ones(size(fmdl.elems,1),1) );
   node_v= calc_all_node_voltages( himg );

   pp= fwd_model_parameters( fmdl );

   elem_v= reshape(node_v( fmdl.elems',:),pp.n_dims+1,[],pp.n_elec);
   elem_v= squeeze(mean(elem_v,1));

   meas_v= pp.N2E * node_v;

   Jbp = zeros(pp.n_elem, sum(fmdl.meas_select) );

   idx= 1;
   % This will only work for stim_patterns with bipolar injection
   for i = 1:pp.n_stim;
     meas_pat_i= fmdl.stimulation(i).meas_pattern;
     elem_vi= elem_v(:,i);
     for j= 1:size(meas_pat_i,1); 
        idx_pl= find(meas_pat_i(j,:)>0);
        idx_mi= find(meas_pat_i(j,:)<0);
        Jbp_idx  = (elem_vi <= meas_v(idx_pl,i)) & ...
                   (elem_vi >  meas_v(idx_mi,i));
        Jbp(:,idx) = - Jbp_idx;
        idx= idx+1;
     end
   end

   if strcmp(type,'naive')
% do nothing
   elseif strcmp(type,'simple_filter')
      xe=mean(reshape(fmdl.nodes(fmdl.elems',1),pp.n_dims+1,pp.n_elem));
      ye=mean(reshape(fmdl.nodes(fmdl.elems',2),pp.n_dims+1,pp.n_elem));
      filt=sqrt(1-xe.^2-ye.^2);
      Jbp= Jbp.*( filt'*ones(1,size(Jbp,2)) );
   else
      error(['dont know what to do with filter type=',type]);
   end

