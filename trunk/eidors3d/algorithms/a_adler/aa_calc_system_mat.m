function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_calc_system_mat.m,v 1.16 2008-05-10 19:12:36 aadler Exp $

p= aa_fwd_parameters( fwd_model );

d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

SSiidx= floor([0:d*e-1]'/d)*d*ones(1,d) + ones(d*e,1)*(1:d) ;
SSjidx= [1:d*e]'*ones(1,d);
SSdata= zeros(d*e,d);
% FIXME: test_2d_resitor gives wrong results (unless dfact = 4*(d-1))
dfact = (d-2)*(d-1); % to match analytic solution 4*dims
for j=1:e
  a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  idx= d*(j-1)+1 : d*j;
  SSdata(idx,1:d)= a(2:d,:)'*a(2:d,:)/abs(det(a));
end %for j=1:ELEMs 

idx= ceil( (1:e*d)/d );
elem_sigma = img.elem_data(idx);
SSdata= SSdata.*( elem_sigma(:) * ones(1,d) );

if 0 % Not complete electrode model
   SS= sparse(SSiidx,SSjidx,SSdata/dfact);
   CC= sparse((1:d*e),p.ELEM(:),ones(d*e,1), d*e, n);
else
   [S2data,S2iidx,S2jidx, C2data,C2iidx,C2jidx] = ...
             compl_elec_mdl(fwd_model,p);
   SS= sparse([SSiidx(:);       S2iidx(:)],...
              [SSjidx(:);       S2jidx(:)],...
              [SSdata(:)/dfact; S2data(:)]);
   
   CC= sparse([(1:d*e)';    C2iidx(:)], ...
              [p.ELEM(:);   C2jidx(:)], ...
              [ones(d*e,1); C2data(:)]);
end

s_mat.E= CC'* SS * CC;


% Add parts for complete electrode model
function [SSdata,SSiidx,SSjidx, CCdata,CCiidx,CCjidx] = ...
             compl_elec_mdl(fwd_model,pp);
   d0= pp.n_dims;
   SSdata= zeros(0,d0);
   SSd_block= ( ones(d0) + eye(d0) )/6/(d0-1); % 6 in 2D, 12 in 3D 
   SSiidx= zeros(0,d0);
   SSjidx= zeros(0,d0);
   SSi_block= ones(d0,1)*(1:d0);
   CCdata= zeros(0,d0);
   CCd_block= ( ones(d0) + eye(d0) )/6/(d0-1); % 6 in 2D, 12 in 3D 
   CCiidx= zeros(0,d0);
   CCjidx= zeros(0,d0);
  
   bdy = double( fwd_model.boundary ); % double if for matlab bugs
   sidx= (d0+1)*pp.n_elem;
   for i= 1:length(fwd_model.electrode);
      zc=  fwd_model.electrode(i).z_contact;
      ffb = find_bdy_idx( bdy, fwd_model.electrode(i).nodes);

      for ff= ffb(:)'
         bdy_nds= bdy(ffb,:);
         bdy_pts= fwd_model.nodes(bdy_nds,:);
         area= tria_area( bdy_pts ); 

         SSdata= [SSdata; SSd_block * area/zc];
         SSiidx= [SSiidx; SSi_block' + sidx];
         SSjidx= [SSjidx; SSi_block  + sidx];

         CCiidx= [CCiidx; SSi_block(1:2,:) + sidx];
         CCjidx= [CCjidx; bdy_nds ; (pp.n_node+i)*ones(1,d0)];
         CCdata= [CCdata; [1;-1]*ones(1,d0)];

         sidx = sidx + d0;
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
