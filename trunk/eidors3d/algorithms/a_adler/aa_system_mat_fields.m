function FC= aa_system_mat_fields( fwd_model )
% AA_SYSTEM_MAT_FIELDS: fields (elem to nodes) fraction of system mat
% FC= aa_system_mat_fields( fwd_model )
% input: 
%   fwd_model = forward model
% output:
%   FC:        s_mat= C' * S * conduct * C = FC' * conduct * FC;

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_system_mat_fields.m,v 1.1 2008-05-11 00:20:55 aadler Exp $

   p= aa_fwd_parameters( fwd_model );
   d0= p.n_dims+0;
   d1= p.n_dims+1;
   e= p.n_elem;
   n= p.n_node;

   FFiidx= floor([0:d0*e-1]'/d0)*d1*ones(1,d1) + ones(d0*e,1)*(1:d1);
   FFjidx= [1:d0*e]'*ones(1,d1);
   FFdata= zeros(d0*e,d1);
   dfact = (d0-1)*d0;
   for j=1:e
     a=  inv([ ones(d1,1), p.NODE( :, p.ELEM(:,j) )' ]);
     idx= d0*(j-1)+1 : d0*j;
     FFdata(idx,1:d1)= a(2:d1,:)/ sqrt(dfact*abs(det(a)));
   end %for j=1:ELEMs 

% FFjidx and FFiidx are inverse labelled
if 0 % Not complete electrode model
   FF= sparse(FFjidx,FFiidx,FFdata);
   CC= sparse((1:d1*e),p.ELEM(:),ones(d1*e,1), d1*e, n);
else
   [F2data,F2iidx,F2jidx, C2data,C2iidx,C2jidx] = ...
             compl_elec_mdl(fwd_model,p);
   FF= sparse([FFjidx(:); F2iidx(:)],...
              [FFiidx(:); F2jidx(:)],...
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
  
   bdy = double( fwd_model.boundary ); % double if for matlab bugs
   sidx= d0*pp.n_elem;
   cidx= (d0+1)*pp.n_elem;
   for i= 1:length(fwd_model.electrode);
      zc=  fwd_model.electrode(i).z_contact;
      ffb = find_bdy_idx( bdy, fwd_model.electrode(i).nodes);

      for ff= ffb(:)'
         bdy_nds= bdy(ff,:);
         bdy_pts= fwd_model.nodes(bdy_nds,:);
         area= tria_area( bdy_pts ); 

         FFdata= [FFdata; FFd_block * sqrt(area/zc)];
         FFiidx= [FFiidx; FFi_block' + sidx];
         FFjidx= [FFjidx; FFi_block  + cidx];

         CCiidx= [CCiidx; FFi_block(1:2,:) + cidx];
         CCjidx= [CCjidx; bdy_nds ; (pp.n_node+i)*ones(1,d0)];
         CCdata= [CCdata; [1;-1]*ones(1,d0)];
         sidx = sidx + d0;
         cidx = cidx + d0;
      end
      
   end

function ffb = find_bdy_idx( bdy, elec_nodes);
   bdy_els = zeros(size(bdy,1),1);
   for nd= unique(elec_nodes);
      bdy_els = bdy_els + any(bdy==nd,2);
   end
   ffb = find(bdy_els == size(bdy,2));

% bdy points is [x1,y1,z1;x2,y2,z2; etc]
function area= tria_area( bdy_pts ); 
   vectors= diff(bdy_pts); 
   if size(vectors,1)==2
      vectors= cross(vectors(1,:),vectors(2,:));
   end
   d= size(bdy_pts,1);
   area= sqrt( sum(vectors.^2) )/( d-1 );
