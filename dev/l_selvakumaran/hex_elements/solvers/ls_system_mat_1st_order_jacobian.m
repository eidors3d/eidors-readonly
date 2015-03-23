function [BB,CC]= ls_system_mat_1st_order_jacobian( fwd_model, img,elem_data)
% SYSTEM_MAT_1ST_ORDER: SS= system_mat_1st_order( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code that uses quad/hex elements
% img.fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1;
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on');
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling SYSTEM_MAT_1ST_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

% check parametrization
if isfield(img,'current_params') && ~isempty(img.current_params) ... 
        && ~strcmp(img.current_params,'conductivity');
    error('system_mat_1st_order does not work for %s',img.current_params);
end
fwd_model.gnd_node=1;
p.nodes= fwd_model.nodes';
p.n_dims=size(p.nodes,1);
p.n_nodes=size(p.nodes,2);
p.elements=fwd_model.elems';
p.n_elem=size(p.elements,2);
p.n_elec=size(fwd_model.electrode,2);
p.boundary=ls_find_boundary(p.elements');
e=p.n_elem;
d0=p.n_dims;
d1=size(p.elements,1);
FFjidx=floor([0:d0*e-1]'/d0)*d1*ones(1,d1)+ones(d0*e,1)*(1:d1);
FFiidx=[1:d0*e]'*ones(1,d1);
BB1=sparse(d1*e,d1*e);
for i=1:d1;
    
   if d0==2;
     FFdata_i=zeros(d0*e,d1);
     P=[ 0.5773 0.5773;...
         -0.5773 0.5773;...
         -0.5773 -0.5773;...
          0.5773 -0.5773];
       for j=1:e;
            p.elements(:,j);
            a=p.nodes(:,p.elements(:,j));
            dN=(1/4)*[1+P(i,2) -1-P(i,2) -1+P(i,2) 1-P(i,2);...
                      1+P(i,1) 1-P(i,1) -1+P(i,1) -1-P(i,1)];
            J=dN*a';
            B=sqrt(det(J))*(inv(J)*dN);
            idx= d0*(j-1)+1 : d0*j;
            FFdata_i(idx,1:d1)= B(:,1:d1);
        end 
    elseif d0==3;
        P=[-0.5773 -0.5773 0.5773;...
            0.5773 -0.5773 0.5773;...
            0.5773 0.5773 0.5773;...
            -0.5773 0.5773 0.5773;...
            -0.5773 -0.5773 -0.5773;...
            0.5773 -0.5773 -0.5773;...
            0.5773 0.5773 -0.5773;...
            -0.5773 0.5773 -0.5773];
       for j=1:e;
            p.elements(:,j);
            a=p.nodes(:,p.elements(:,j));
            dN=(1/8)*[-(1-P(i,2))*(1-P(i,3)) (1-P(i,2))*(1-P(i,3)) (1+P(i,2))*(1-P(i,3)) -(1+P(i,2))*(1-P(i,3))...
                -(1-P(i,2))*(1+P(i,3)) (1-P(i,2))*(1+P(i,3)) (1+P(i,2))*(1+P(i,3)) -(1+P(i,2))*(1+P(i,3));...
                -(1-P(i,1))*(1-P(i,3)) -(1+P(i,1))*(1-P(i,3)) (1+P(i,1))*(1-P(i,3)) (1-P(i,1))*(1-P(i,3))...
                -(1-P(i,1))*(1+P(i,3)) -(1+P(i,1))*(1+P(i,3)) (1+P(i,1))*(1+P(i,3)) (1-P(i,1))*(1+P(i,3));...
                -(1-P(i,1))*(1-P(i,2)) -(1+P(i,1))*(1-P(i,2)) -(1+P(i,1))*(1+P(i,2)) -(1-P(i,1))*(1+P(i,2))...
                (1-P(i,1))*(1-P(i,2)) (1+P(i,1))*(1-P(i,2)) (1+P(i,1))*(1+P(i,2)) (1-P(i,1))*(1+P(i,2))];
            J=dN*a';
            B=sqrt(det(J))*(inv(J)*dN);
            idx= d0*(j-1)+1 : d0*j;
            FFdata_i(idx,1:d1)= B(:,1:d1);
        end
   end
   FF= sparse(FFiidx(:),...
              FFjidx(:),...
              FFdata_i(:));
 
   BB0=system_matrix_construct(fwd_model,img,FF,elem_data); 
   BB1=BB1+BB0;
                      
end
[F2data,F2iidx,F2jidx, C2data,C2iidx,C2jidx] = ...
             compl_elec_mdl(fwd_model,p);
             
         F2= sparse(F2iidx(:),...
              F2jidx(:),...
              F2data(:));
      
         BB2=F2'*F2;
         [row1,col1,v1] = find(BB1);
         [row2,col2,v2] = find(BB2);
         BB=sparse([row1;row2],[col1;col2],[v1;v2]);
         
         CC= sparse([(1:d1*e)';    C2iidx(:)], ...
              [p.elements(:);   C2jidx(:)], ...
              [ones(d1*e,1); C2data(:)]);



function [FFdata,FFiidx,FFjidx, CCdata,CCiidx,CCjidx] = ...
             compl_elec_mdl(fwd_model,pp);
   d= pp.n_dims;
   if d==2;
    d0=2;
    cidx= 4*pp.n_elem;
    FFd_block= sqrtm( ( ones(d0) + eye(d0) )/6);
   elseif d==3;
    d0=4;   
    cidx= 8*pp.n_elem;
    FFd_block= sqrtm( [1 0.5 0.25 0.5;0.5 1 0.5 0.25;0.25 0.5 1 0.5; 0.5 0.25 0.5 1]/9);
   end
   FFdata= zeros(0,d0);   
   FFiidx= zeros(0,d0);
   FFjidx= zeros(0,d0);
   FFi_block= ones(d0,1)*(1:d0);
   CCdata= zeros(0,d0);
   CCiidx= zeros(0,d0);
   CCjidx= zeros(0,d0);
   sidx= d0*pp.n_elem;
   
   for i= 1:pp.n_elec;
      eleci = fwd_model.electrode(i);
      zc=  eleci.z_contact;
      [bdy_idx, bdy_area] = ls_find_electrode_bdy( ...
          pp.boundary, fwd_model.nodes, eleci.nodes );
      for j= 1:length(bdy_idx);
         bdy_nds= pp.boundary(bdy_idx(j),:);

         FFdata= [FFdata; FFd_block * sqrt(bdy_area(j)/zc)];
         FFiidx= [FFiidx; FFi_block' + sidx];
         FFjidx= [FFjidx; FFi_block  + cidx];

         CCiidx= [CCiidx; FFi_block(1:2,:) + cidx];
         CCjidx= [CCjidx; bdy_nds ; (pp.n_nodes+i)*ones(1,d0)];
         CCdata= [CCdata; [1;-1]*ones(1,d0)];
         sidx = sidx + d0;
         cidx = cidx + d0;
      end
      
   end



function [BB]=system_matrix_construct(fwd_model,img,FF,elem_data);
lFC= size(FF,1);
ED = ls_elem_dim(fwd_model);
lNE= ED*num_elems(fwd_model);

elem_data = check_elem_data(fwd_model, img,elem_data);
if size(elem_data,3) == 1;
% Scalar conductivity == isotropic
   elem_sigma = kron( elem_data, ones(ED,1) );
   ES= spdiags(elem_sigma,0,lFC,lFC);
else
    switch ls_elem_dim(fwd_model);
     case 2;
       idx = 1:2:lNE;
       ES= sparse([idx,idx+1,idx,idx+1]', ...
                  [idx,idx,idx+1,idx+1]', elem_data(:), lFC,lFC);     
     case 3;
       idx = 1:3:lNE;
       ES= sparse([idx,idx+1,idx+2,idx,idx+1,idx+2,idx,idx+1,idx+2]', ...
                  [idx,idx,idx,idx+1,idx+1,idx+1,idx+2,idx+2,idx+2]', ...
                   elem_data(:), lFC,lFC);
     otherwise; 
       error('%d D anisotropic elements not implemented', elem_dim(fwd_model));
    end  
end

BB= FF' * ES * FF;



function elem_data = check_elem_data(fwd_model, img,elem_data);
  % elem_data = img.elem_data; 
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('system_mat_1st_order: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(fwd_model, 'coarse2fine');
     c2f = fwd_model.coarse2fine;
     sz_c2f = size(c2f);
     switch sz_elem_data(1);
       case sz_c2f(1); % Ok     
       case sz_c2f(2); 
           
      if size(elem_data,3)==1;
%     fwd_model data is provided on coarse mesh
      elem_data = c2f * elem_data; 
      else
          switch size(elem_data,3);    
          case 2;
              for i=1:2;
                  for j=1:2;
                      C=c2f*elem_data(:,1,i,j);
                      conduct1(:,1,i,j)=C;
                  end
              end
              elem_data=conduct1;
          case 3;
             for i=1:3;
                  for j=1:3;
                      C=c2f*elem_data(:,1,i,j);
                      conduct1(:,1,i,j)=C;
                  end
             end
              elem_data=conduct1;
          end
       end  
          
           
           if isfield(fwd_model, 'background');
              elem_data = elem_data + fwd_model.background;
          end

       otherwise; error(['system_mat_1st_order: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(fwd_model);
       error(['system_mat_1st_order: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end

