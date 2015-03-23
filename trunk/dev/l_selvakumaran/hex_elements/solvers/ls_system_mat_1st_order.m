function s_mat= ls_system_mat_1st_order( fwd_model, img)
% LS_SYSTEM_MAT_1ST_ORDER: SS= system_mat_1st_order( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code which uses quad elements
% img.fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1;
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
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

%fwd_model.electrode=[];
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
FFjidx=floor((0:d0*e-1)'/d0)*d1*ones(1,d1)+ones(d0*e,1)*(1:d1);
FFiidx=(1:d0*e)'*ones(1,d1);
FC=[];
E=sparse(p.n_nodes,p.n_nodes);
s_mat.E=sparse(p.n_nodes+p.n_elec,p.n_nodes+p.n_elec);

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
            B=sqrt(det(J))*(J\dN);
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
            B=sqrt(det(J))*(J\dN);
            idx= d0*(j-1)+1 : d0*j;
            FFdata_i(idx,1:d1)= B(:,1:d1);
        end
   end
   FF= sparse(FFiidx(:),...
              FFjidx(:),...
              FFdata_i(:));
   
   CC= sparse((1:d1*e)', ...
              p.elements(:), ...
              ones(d1*e,1));

   FC= FF*CC;
   s_mat_i=system_matrix_construct(fwd_model,img,FC);
   E=E+s_mat_i.E;
                      
end
[full_E] = compl_elec_mdl(fwd_model,p,E);
s_mat.E=full_E;



function [Ef] = compl_elec_mdl(fwd_model,pp,E)
   vr=pp.n_nodes;
   er=pp.n_elec;
   E=full(E);
   Ef=spalloc(vr+er,vr+er,er*vr);
   Ef(1:vr,1:vr)=E;
   
for i=1:pp.n_elec;
       td=0;
       
       eleci = fwd_model.electrode(i);
       zc=  eleci.z_contact;
       [bdy_idx, bdy_area] = ls_find_electrode_bdy( ...
          pp.boundary, pp.nodes', eleci.nodes );
     if pp.n_dims==2; 
       for j= 1:length(bdy_idx);
            bdy_nodes   = pp.boundary(bdy_idx(j),:);
            m=bdy_nodes(1,1);
            n=bdy_nodes(1,2);

           td=td+bdy_area(j);

           Ef(m,vr+i) = Ef(m,vr+i) - bdy_area(j)/(zc*2) ; % Kv -> Ec  -> Vertical bar
           Ef(n,vr+i) = Ef(n,vr+i) - bdy_area(j)/(zc*2) ; % Kv -> Ec

           Ef(vr+i,m) = Ef(vr+i,m) - bdy_area(j)/(zc*2) ; % Kv' -> Ec' -> Horizontal bar
           Ef(vr+i,n) = Ef(vr+i,n) - bdy_area(j)/(zc*2) ; % Kv' -> Ec'


           Ef(m,m) = Ef(m,m) + bdy_area(j)/(zc*3); % Kz -> E -> Main bar
           Ef(n,n) = Ef(n,n) + bdy_area(j)/(zc*3); % Kz -> E
           Ef(m,n) = Ef(m,n) + bdy_area(j)/(zc*6); % Kz -> E
           Ef(n,m) = Ef(n,m) + bdy_area(j)/(zc*6); % Kz -> E

       end
       Ef(vr+i,vr+i) = Ef(vr+i,vr+i) + td/(zc);% dealing with this electrode
     elseif pp.n_dims==3;
        for j= 1:length(bdy_idx);
            bdy_nodes   = pp.boundary(bdy_idx(j),:);
            m=bdy_nodes(1,1);
            n=bdy_nodes(1,2);
            o=bdy_nodes(1,3);
            p=bdy_nodes(1,4);

           td=td+bdy_area(j);

           Ef(m,vr+i) = Ef(m,vr+i) - bdy_area(j)/(4*zc) ; % Kv -> Ec  -> Vertical bar
           Ef(n,vr+i) = Ef(n,vr+i) - bdy_area(j)/(4*zc) ; % Kv -> Ec
           Ef(o,vr+i) = Ef(o,vr+i) - bdy_area(j)/(4*zc) ;
           Ef(p,vr+i) = Ef(p,vr+i) - bdy_area(j)/(4*zc) ;
           
           
           Ef(vr+i,m) = Ef(vr+i,m) - bdy_area(j)/(4*zc) ; % Kv' -> Ec' -> Horizontal bar
           Ef(vr+i,n) = Ef(vr+i,n) - bdy_area(j)/(4*zc) ; % Kv' -> Ec'
           Ef(vr+i,o) = Ef(vr+i,o) - bdy_area(j)/(4*zc) ;
           Ef(vr+i,p) = Ef(vr+i,p) - bdy_area(j)/(4*zc) ;

           Ef(m,m) = Ef(m,m) + bdy_area(j)/(zc*9); % Kz -> E -> Main bar
           Ef(n,n) = Ef(n,n) + bdy_area(j)/(zc*9); 
           Ef(o,o) = Ef(o,o) + bdy_area(j)/(zc*9);
           Ef(p,p) = Ef(p,p) + bdy_area(j)/(zc*9);% Kz -> E
           Ef(m,n) = Ef(m,n) + bdy_area(j)/(2*zc*9); % Kz -> E
           Ef(n,m) = Ef(n,m) + bdy_area(j)/(2*zc*9); % Kz -> E
           Ef(m,o) = Ef(m,o) + bdy_area(j)/(4*zc*9);
           Ef(o,m) = Ef(o,m) + bdy_area(j)/(4*zc*9);
           Ef(m,p) = Ef(m,p) + bdy_area(j)/(2*zc*9);
           Ef(p,m) = Ef(p,m) + bdy_area(j)/(2*zc*9);
           Ef(n,o) = Ef(n,o) + bdy_area(j)/(2*zc*9);
           Ef(o,n) = Ef(o,n) + bdy_area(j)/(2*zc*9);
           Ef(n,p) = Ef(n,p) + bdy_area(j)/(4*zc*9);
           Ef(p,n) = Ef(p,n) + bdy_area(j)/(4*zc*9);
           Ef(o,p) = Ef(o,p) + bdy_area(j)/(2*zc*9);
           Ef(p,o) = Ef(p,o) + bdy_area(j)/(2*zc*9);
        end% dealing with this electrode
        Ef(vr+i,vr+i) = Ef(vr+i,vr+i) + td/(zc);
     end
end


function [s_mat_i]=system_matrix_construct(fwd_model,img,FC);
lFC= size(FC,1);
ED = ls_elem_dim(fwd_model);
lNE= ED*num_elems(fwd_model);

elem_data = check_elem_data(fwd_model, img);
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

s_mat_i.E= FC' * ES * FC;



function elem_data = check_elem_data(fwd_model, img)
   elem_data = img.elem_data; 
   sz_elem_data = size(elem_data);
   size(elem_data,3);
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
      
           
          if isfield(fwd_model, 'background')
              elem_data = elem_data + fwd_model.background;
          end

       otherwise; error(['system_mat_1st_order: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(fwd_model)
       error(['system_mat_1st_order: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end

