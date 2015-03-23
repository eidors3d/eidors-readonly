function [J]= ls_jacobian_adjoint( img)
% AA_CALC_JACOBIAN: J= jacobian_adjoint_simplecode( fwd_model, img)
% Calculate Jacobian Matrix for EIT Alg of Adler & Guardo 1996 for quad/hex
% elements
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id$

fwd_model = img.fwd_model;
pp= ls_fwd_model_parameters( fwd_model );
s_mat= ls_system_mat_1st_order( fwd_model, img );
d= size(pp.ELEM,1);
if ~isfield(fwd_model,'coarse2fine');
 e= pp.n_elem;
else 
 e= size(fwd_model.coarse2fine,2);
end
n= pp.n_node;

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
DE= zeros(pp.n_elec, pp.n_stim, e);

idx= 1:size(s_mat.E,1);
idx( fwd_model.gnd_node ) = [];
sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat.E(idx,idx) \ pp.QQ( idx,: );

zi2E(:,idx)= pp.N2E(:,idx)/ s_mat.E(idx,idx);

if size(img.elem_data,3)==1;
aniso=1;
jacob.elem_data=img.elem_data;
else
    switch size(pp.NODE,1);
        case 2;
          if ~all(img.elem_data(:,1,1,2)) &&  ~all(img.elem_data(:,1,2,1));
              aniso=2;
              jacob.elem_data=zeros(e,aniso,2,2);
              for i=1:aniso;
                jacob.elem_data(:,i,i,i)=ones(e,1);
              end
          else 
              aniso=3;
              jacob.elem_data=zeros(e,aniso,2,2);
              for i=1:2;
              jacob.elem_data(:,i,i,i)=ones(e,1);
              end
              jacob.elem_data(:,3,1,2)=ones(e,1);
              jacob.elem_data(:,3,2,1)=ones(e,1);
          end
        case 3;
            if ~all(img.elem_data(:,1,1,2)) &&  ~all(img.elem_data(:,1,2,1)) && ~all(img.elem_data(:,1,3,1)) &&...
                ~all(img.elem_data(:,1,1,3)) && ~all(img.elem_data(:,1,2,3)) && ~all(img.elem_data(:,1,3,2));
              aniso=3;
              jacob.elem_data=zeros(e,aniso,3,3);
              for i=1:aniso;
                jacob.elem_data(:,i,i,i)=ones(e,1);
              end
          else 
              aniso=6;
              jacob.elem_data=zeros(e,aniso,3,3);
              for i=1:3;
                jacob.elem_data(:,i,i,i)=ones(e,1);
              end
              jacob.elem_data(:,4,1,2)=ones(e,1);
              jacob.elem_data(:,5,2,3)=ones(e,1);
              jacob.elem_data(:,6,3,1)=ones(e,1);
              jacob.elem_data(:,4,2,1)=ones(e,1);
              jacob.elem_data(:,5,3,2)=ones(e,1);
              jacob.elem_data(:,6,1,3)=ones(e,1);
            end
    end
end

%----------------------------------------------------------------------------
DE_an=zeros(pp.n_elec, pp.n_stim, 0);
for i=1:aniso;
[dSS_dEj,CC]= ls_system_mat_1st_order_jacobian( fwd_model, img,jacob.elem_data(:,i,:,:));
if ~isfield(fwd_model,'coarse2fine');
for k= 1:e;
    idx= d*(k-1)+1 : d*k;
    CC_idx = CC(idx,:);
    dSS_dEj1=dSS_dEj(idx,idx);
    dq= zi2E * CC_idx' * dSS_dEj1 * CC_idx * sv;
    DE(:,:,k)= dq;
end
else
    c2f=fwd_model.coarse2fine;
    DE= zeros(pp.n_elec, pp.n_stim, size(c2f,2) );
      for k= 1:e;
          ff = find( c2f(:,k) );
          lff= length(ff)*(d);
          ff1= ones(d,1) * ff(:)';
          ffd= (d)*ff1 + (-(d-1):0)'*ones(1,length(ff));
          dDD_dEj = spdiags(c2f(ff1,k), 0, lff, lff);
          CC_idx=CC(ffd,:);
          dSS_dEj1=dSS_dEj(ffd,ffd); 
          dq= zi2E * CC_idx' * dSS_dEj1 * dDD_dEj * CC_idx * sv;
          DE(:,:,k)= dq;
      end
end
    DE_an=cat(3,DE_an,DE);
end


J = zeros( pp.n_meas, e*aniso );
idx=0;
for j= 1:pp.n_stim;
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   DEj = reshape( DE_an(:,j,:), pp.n_elec, e*aniso );
   J( idx+(1:n_meas),: ) = meas_pat*DEj;
   idx= idx+ n_meas;
end

% calculate normalized Jacobian
if pp.normalize;
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,e));
   
end

% Calculated Jacobian is inverted
J = -J;
