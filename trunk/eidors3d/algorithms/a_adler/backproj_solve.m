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
% $Id: backproj_solve.m,v 1.3 2007-11-28 15:07:39 aadler Exp $

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
    Jbp = calc_backprojection_mask( fwd_model, type );

    one_step_inv= Jbp;

    eidors_obj('set-cache', inv_model, 'backproj_solve', one_step_inv);
    eidors_msg('backproj_solve: setting cached value', 2);
end

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

sol = one_step_inv * dv;

% create a data structure to return
img.name= 'solved by backproj_solve';
img.elem_data = sol;
img.fwd_model= fwd_model;

function Jbp = calc_backprojection_mask( fmdl , type);
   % create homog image
   himg= eidors_obj('image','', 'fwd_model',fmdl, ...
                    'elem_data',ones(size(fmdl.elems,1),1) );
   node_v= calc_all_node_voltages( himg );

   pp= aa_fwd_parameters( fmdl );

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
        idx_pl= find(meas_pat_i(j,:)<0);
        idx_mi= find(meas_pat_i(j,:)>0);
        Jbp(:,idx) = (elem_vi <= meas_v(idx_pl,i)) & ...
                     (elem_vi >  meas_v(idx_mi,i));
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

