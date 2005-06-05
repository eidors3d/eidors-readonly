function J= aa_calc_jacobian( fwd_model, img)
% AA_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for EIT Alg of Adler & Guardo 1996
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc
% $Id: aa_calc_jacobian.m,v 1.1 2005-06-05 13:43:31 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

J= zeros( pp.n_meas, e );


idx= 1:pp.n_node;
idx( fwd_model.gnd_node ) = [];
sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat(idx,idx) \ pp.QQ;


-----------------
% measured voltages from v
vv = zeros( pp.n_meas, 1 );
idx=0;
for i=1:pp.n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,2);
   vv( idx+(1:n_meas) ) = meas_pat'*v_els(:,i);
   idx= idx+ n_meas;
end


  i=1:n; i(volt(1,:))=[];
  DVV= zeros(sum(ELS),c);
  z= CC'*SS*CC;
  zinv=zeros(n);
  zinv(i,i)= inv(z(i,i));
  end
  dzi= zinv(MES(1,:),:)- zinv(MES(2,:),:);
  vvh= dzi*QQ;
  sv= zinv*QQ;
  for j=1:c
    idx= find(ones(d,1)*(CPTR==j));
    dq= dzi*(  CC(idx,:)'*SS(idx,idx)*CC(idx,:)  )*sv;
    DVV(:,j)= dq(ELS)./vvh(ELS);
  end
