function FC= aa_system_mat_fields( fwd_model )
% AA_SYSTEM_MAT_FIELDS: fields (elem to nodes) fraction of system mat
% FC= aa_system_mat_fields( fwd_model )
% input: 
%   fwd_model = forward model
% output:
%   FC:        s_mat= C' * S * conduct * C = FC' * conduct * FC;

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

cache_obj = mk_cache_obj(fwd_model);
FC = eidors_obj('get-cache', cache_obj, 'aa_system_mat_fields');
if ~isempty(FC)
   eidors_msg('aa_system_mat_fields: using cached value', 4);
   return
end

FC= calc_system_mat_fields( fwd_model );

eidors_cache('boost_priority',1); % Moderate Priority boost
eidors_obj('set-cache', cache_obj, 'aa_system_mat_fields', FC);
eidors_msg('aa_system_mat_fields: setting cached value', 4);
eidors_cache('boost_priority',-1);

% only cache stuff which is really relevant here
function cache_obj = mk_cache_obj(fwd_model);
   cache_obj.elems       = fwd_model.elems;
   cache_obj.nodes       = fwd_model.nodes;
   try
   cache_obj.electrode   = fwd_model.electrode; % if we have it
   end
   cache_obj.type        = 'fwd_model';
   cache_obj.name        = ''; % it has to have one

function FC= calc_system_mat_fields( fwd_model );
   p= aa_fwd_parameters( fwd_model );
   d0= p.n_dims+0;
   d1= p.n_dims+1;
   e= p.n_elem;
   n= p.n_node;

   FFjidx= floor([0:d0*e-1]'/d0)*d1*ones(1,d1) + ones(d0*e,1)*(1:d1);
   FFiidx= [1:d0*e]'*ones(1,d1);
   FFdata= zeros(d0*e,d1);
   dfact = (d0-1)*d0;
   for j=1:e
     a=  inv([ ones(d1,1), p.NODE( :, p.ELEM(:,j) )' ]);
     idx= d0*(j-1)+1 : d0*j;
     FFdata(idx,1:d1)= a(2:d1,:)/ sqrt(dfact*abs(det(a)));
   end %for j=1:ELEMs 

if 0 % Not complete electrode model
   FF= sparse(FFiidx,FFjidx,FFdata);
   CC= sparse((1:d1*e),p.ELEM(:),ones(d1*e,1), d1*e, n);
else
   [F2data,F2iidx,F2jidx, C2data,C2iidx,C2jidx] = ...
             compl_elec_mdl(fwd_model,p);
   FF= sparse([FFiidx(:); F2iidx(:)],...
              [FFjidx(:); F2jidx(:)],...
              [FFdata(:); F2data(:)]);
   
   CC= sparse([(1:d1*e)';    C2iidx(:)], ...
              [p.ELEM(:);   C2jidx(:)], ...
              [ones(d1*e,1); C2data(:)]);
end

FC= FF*CC;

% Add parts for complete electrode model
function [FFdata,FFiidx,FFjidx, CCdata,CCiidx,CCjidx] = ...
             compl_elec_mdl(fwd_model,pp);
   d0= pp.n_dims;
   FFdata= zeros(0,d0);
   FFd_block= sqrtm( ( ones(d0) + eye(d0) )/6/(d0-1) ); % 6 in 2D, 12 in 3D 
   FFiidx= zeros(0,d0);
   FFjidx= zeros(0,d0);
   FFi_block= ones(d0,1)*(1:d0);
   CCdata= zeros(0,d0);
   CCiidx= zeros(0,d0);
   CCjidx= zeros(0,d0);
  
   sidx= d0*pp.n_elem;
   cidx= (d0+1)*pp.n_elem;
   for i= 1:pp.n_elec
      zc=  fwd_model.electrode(i).z_contact;
%     ffb = find_bdy_idx( bdy, fwd_model.electrode(i).nodes);
      [bdy_idx, bdy_area] = find_electrode_bdy( bdy, fwd_model.nodes, ...
                               fwd_model.electrode(i).nodes);

      for j= 1:length(bdy_idx);
         bdy_nds= pp.boundary(bdy_idx(j),:);

         FFdata= [FFdata; FFd_block * sqrt(bdy_area(j)/zc)];
         FFiidx= [FFiidx; FFi_block' + sidx];
         FFjidx= [FFjidx; FFi_block  + cidx];

         CCiidx= [CCiidx; FFi_block(1:2,:) + cidx];
         CCjidx= [CCjidx; bdy_nds ; (pp.n_node+i)*ones(1,d0)];
         CCdata= [CCdata; [1;-1]*ones(1,d0)];
         sidx = sidx + d0;
         cidx = cidx + d0;
      end
      
   end
