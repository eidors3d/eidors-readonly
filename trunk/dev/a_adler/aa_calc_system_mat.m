function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id$

p= fwd_model_parameters( fwd_model );

d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;
SS = sparse(d*e,d*e);

dfact= factorial(d-1);
for j=1:e
  a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  area = 1/dfact/abs(det(a));
  idx= d*(j-1)+1 : d*j;
  SS(idx,idx)= area*a(2:d,:)'*a(2:d,:);
end %for j=1:ELEMs 

CC= sparse((1:d*e),p.ELEM(:),ones(d*e,1), d*e, n);

idx= 1:e*d;
if 1
elem_sigma = sparse(idx,idx, img.elem_data(ceil(idx/d)) );
s_mat.E= CC'* SS * elem_sigma * CC;
end

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
      eleci = fwd_model.electrode(i);
      zc=  eleci.z_contact;
%     ffb = find_bdy_idx( bdy, fwd_model.electrode(i).nodes);
      [bdy_idx, bdy_area] = find_electrode_bdy( ...
          pp.boundary, fwd_model.nodes, eleci.nodes );

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

