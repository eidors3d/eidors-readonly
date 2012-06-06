function HP= choose_noise_figure( inv_model );
% CHOOSE_NOISE_FIGURE: choose hyperparameter based on NF calculation
% HP= choose_noise_figure( inv_model );
% inv_model  => inverse model struct
%
% In order to use this function, it is necessary to specify
% inv_model.hyperparameter. has the following fields
% hpara.func         = @choose_noise_figure;
% hpara.noise_figure = NF Value requested
% hpara.tgt_elems    = vector of element numbers of contrast in centre
%
% The NF parameter is from the definition in Adler & Guardo (1996).
%
% SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
% SNR_x = sumsq(A*x0) / var(A*x) = sum((A*x0).^2) / trace(ABRnB'A')
%   where Rn = noise covariance and A_ii = area of element i
% NF = SNR_z / SNR_x

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

reqNF= inv_model.hyperparameter.noise_figure;
try % remove the value field
   inv_model.hyperparameter = rmfield(inv_model.hyperparameter,'value');
end

HP = eidors_obj('get-cache', inv_model, 'HP_for_NF');
if ~isempty(HP)
   eidors_msg('choose_noise_figure: using cached value', 2);
   return
end

HP = HP_for_NF_search(reqNF,inv_model);

   
eidors_cache('boost_priority',3);
eidors_obj('set-cache', inv_model, 'HP_for_NF', HP);
eidors_msg('choose_noise_figure: setting cached value', 2);
eidors_cache('boost_priority',-3);


function [HP,NF,SE] = HP_for_NF_search(dNF,imdl);
   hp= search1(dNF, imdl, 1);

   dx= hp-linspace(-0.7,0.7,5);
   hp= search2(dNF, imdl, dx);

   dx= hp-linspace(-0.2,0.2,5);
   hp= search2(dNF,imdl, dx);

   dx= hp-linspace(-0.1,0.1,5); 
   hp= search2(dNF,imdl, dx);

   dx= hp-linspace(-0.05,0.05,21);
   hp= search2(dNF,imdl, dx);

   HP= 10^-hp;
  
function hp= search1(dNF, imdl, hp)
  [NF,SE]=calc_noise_figure( imdl, 10^(-hp));
   if     NF+3*SE < dNF; dir = 1;
   elseif NF-3*SE > dNF; dir = -1;
   else   dir = 0; end
   while  dir*NF+3*SE < dir*dNF %>
     hp= hp+0.5*dir;
     [NF,SE]=calc_noise_figure( imdl, 10^(-hp));
   end

function hp=search2(dNF, imdl, dx)
   hp = [];
   it = 0;
   nf = zeros(1,length(dx));
   while isempty(hp) && it < 2
       it = it+1;
       %if it > 1 keyboard, end
       for k=1:length(dx)
           nf(k)=nf(k)+calc_noise_figure( imdl, 10^-dx(k), 10 );
       end
       log_nf = log10(nf/it);
       p= polyfit( dx, log_nf-log10(dNF), 3);
       hp = roots(p);
       %eliminate complex roots
       hp = hp(imag(hp)==0);
       %eliminate out of range roots
       hp = hp( hp<max(dx) & hp>min(dx) );  %USE if poly>1
       % pick the root with the smallest derivative
       if numel(hp) >1
           p2 = p.*(numel(p)-1:-1:0); p2(end) = [];
           [jnk, pos] = min(abs(polyval(p2,hp)));
           hp = hp(pos);
       end
   end
   if isempty(hp)
       %fallback
       [jnk,idx] = min(abs(log_nf-log10(dNF)));
       hp = dx(idx);
   end
