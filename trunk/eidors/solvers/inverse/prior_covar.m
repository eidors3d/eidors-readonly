function Reg= prior_covar( inv_model )
% PRIOR_COVAR image prior with distance-based interelement covar
% This is a simplification of prior_exponential_covar.m
% Reg= prior_covar( inv_model )
% Reg        => output regularization 
% inv_model  => inverse model struct
% P_type--prior type
% 1: elements are globally correlated
% 2: elements within/without electrode rings are correlated to elements in same region.
% 3: only elements within electrode rings are correlated.

% (C) 2007, Tao Dai and Andy Adler. Licenced under the GPL Version 2
% $Id$

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end
cache_obj = {inv_model.fwd_model, inv_model.fourD_prior.P_type};
Reg = eidors_obj('get-cache', cache_obj, 'prior_covar');
if ~isempty(Reg)
   eidors_msg('inv_solve_4d_prior: using cached value', 4);
   return;
else
   Reg = prior_covar_calc( cache_obj);
   eidors_obj('set-cache', cache_obj, 'prior_covar', Reg);
   eidors_msg('prior_covar: setting cached value', 4);
end

function Reg = prior_covar_calc( cache_obj );
   ff = cache_obj{1}; 
   P_type = cache_obj{2};
   % get average x,y,z of each element
   nel = size(ff.elems,1);
   eta = .1;%attenuation factor. eta is large when elems're spatially highly correlated

   dist = zeros(nel);

   z1 = ff.nodes(ff.electrode(1).nodes,3);%upper electrode ring
   z2 = ff.nodes(ff.electrode(end).nodes,3);%lower electrode ring

   for dim = 1: size(ff.nodes,2);
       coords = reshape(ff.nodes(ff.elems,dim),[size(ff.elems)]);%coord of four vertex of each elem
       m_coords = mean( coords,2);%calc center coord
       if dim == 3
           temp = double((m_coords<max(z1,z2))&(m_coords>min(z1,z2)));%regions of interests
           H = temp*temp';
           switch P_type %Prior type
               case 1 % all elements are correlated
                   H = eta*(ones(size(H)));
               case 2 % elements of ROI are not correlated with elements of non_ROI
                   temp = double(m_coords>max(z1,z2));
                   H = H + temp*temp';
                   temp = double(m_coords<min(z1,z2));
                   H = H + temp*temp';
                   H = eta*(H+1e-6);
               case 3 % elements are only correlated within ROI, non-ROI elements are highly attenuated
                   H = eta*(H+1e-6);
               otherwise
                   error('no such a 3-D prior type');
           end
       end
       difm = m_coords*ones(1,nel);
       difm = difm - difm';
       dist= dist + difm.^2;%
   end
   dist = sqrt(dist);%elements distance matrix
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   dist = dist/max(dist(:));

   Reg = exp(-dist ./ H);% 3-D elements correlations matrix.

function do_unit_test
   imdl = mk_common_model('b3cr',[16,3]);
   imdl.fourD_prior.P_type = 2;
tic; 
   prior_covar( imdl );
toc;
