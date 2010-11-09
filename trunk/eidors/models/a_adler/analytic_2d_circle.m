function V = analytic_2d_circle(I, params)
% V = analytic_2d_circle(I, [s_h, s_i, b, a, angl])
% Voltage around a 2D circle with properties
%   conductivities: tank(s_h), inclusion(s_i)
%   tank radius = 1
%   inclusion: radius(a), distance(b) at angle(angl)
%
% I = current at a regular sample of nodes on the boundary
%
% Based on eqn 21 in Seagar and Bates (1985).
%
% (C) Andy Adler 2009. License: GPL V2 or V3
% $Id$
%
% USAGE EXAMPLE: compare FEM to analytic model
%
%  imdl= mk_common_model('c2c0',16);
% img = calc_jacobian_bkgnd(imdl);
% img.fwd_solve.get_all_meas = 1;
% img.elem_data(1:256)= 0.1;
% vi= fwd_solve(img);
% 
% [jnk,maxl] = min(vi.volt(:,1));
% vi= vi.volt(maxl:end,1);
% lv= length(vi);
% 
% vol = get_elem_volume(img.fwd_model);
% vt= sum(vol); va= sum(vol(img.elem_data ~=1));
% rad = sqrt(va/vt);
% 
% I =  zeros(lv,1); I(1+[0,lv/16]) = [-25,25]*lv/16;
% vsi= analytic_2d_circle(I, [1, 0.1, 0.0, rad, 0]);
% 
% plot([vi,vsi]);


if ~all(abs(sum(I,1)) <1e-12); error('net I must be 0'); end
[ll,nv] = size(I);

IF= fft(I); IF(1)= []; % Don't include zero term
[SM,ll2] = solve_matrix(ll, params);
VF2= SM*IF(ll2);
VF= zeros(size(I));
VF(1+[ll2, ll-ll2],:)= [VF2;conj(VF2)];
V = ifft(VF);

if norm(imag(V))>1e-10; error('Unexpected Imaginary output - probably a bug'); end
V= real(V);


function [SM,ll2] = solve_matrix(ll, params)
  ll1= ll-1; % IGNORE DC
  ll2 = 1:ceil(ll1/2);

  SM = eidors_obj('get-cache', {ll,params}, 'analytic_2d_circle');
  if ~isempty(SM)
     eidors_msg('analytic_2d_circle: using cached value', 3);
     return
  end
  
   s_h = params(1);
   s_i = params(2);
   b   = params(3);
   a   = params(4);
   alp = params(5);

   mll2= max(ll2);
   K = ll2(:);
   D = spdiags(K,0,mll2,mll2);
   iD = spdiags(1./K,0,mll2,mll2); % inv D

   if s_h == s_i 
      SM= 1/s_h * iD;
      return
   end

   mu = (s_h - s_i) / (s_h + s_i);

   if b==0
      mua2k = mu*a.^(2*K);
      K = K.*(1 - mua2k)./(1 + mua2k);
      iD = spdiags(1./K,0,mll2,mll2); % inv(D)

      SM= 1/s_h * iD;
      return
   end



   T = sparse(mll2,mll2);
   for m=ll2; for n=ll2
     sm= 0;
     for p=1:min(m,n)
%      sm=sm+ nchoosek(m-1,p-1)*nchoosek(n,p)*ab2^p;
       % stop calculating nchoosek when too small, because 
       dsm = nchoosek(m-1,p-1)*nchoosek(n,p)*a^(2*p)*b^(m+n-2*p);
%       disp([m,n,p,log10(dsm)])
       if dsm < 1e-11; break; end  
       sm=sm + dsm;
     end
%    T(i,j) = mu*(sm/m)* b^(m+n) * exp(1j*(m-n)*alp);
     T(m,n) = mu*(sm/m)* exp(1j*(m-n)*alp);
   end; end

   I= speye(mll2);
   SM = (I-T*D)\(I+T*D)*iD;

   eidors_obj('set-cache', {ll,params}, 'analytic_2d_circle', SM);
   eidors_msg('analytic_2d_circle: setting cached value', 3);