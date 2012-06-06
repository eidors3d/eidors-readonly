function Reg= tutorial210_cheat_laplace( inv_model )
% Reg= cheat_laplace( inv_model )
% Reg        => output regularization term
% Parameters:
%   elems    = inv_model.tutorial210_cheat_laplace.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.tutorial210_cheat_laplace.cheat_weight;
%            => new weight to set elements to

pp= fwd_model_parameters( inv_model.fwd_model );

ROI = zeros(1,pp.n_elem);
ROI( inv_model.tutorial210_cheat_laplace.cheat_elements ) = 1;

Iidx= [];
Jidx= [];
Vidx= [];
for ii=1:pp.n_elem
  el_adj = find_adjoin( ii, pp.ELEM );
  for jj=el_adj(:)'
      if (ROI(ii) + ROI(jj)) == 1 %one only
         fac= inv_model.tutorial210_cheat_laplace.cheat_weight *.5;
      else 
         fac = .5;
      end
      Iidx= [Iidx,      ii, ii, jj, jj];
      Jidx= [Jidx,      ii, jj, ii, jj];
      Vidx= [Vidx, fac*([1, -1, -1,  1]) ];
  end
end
Reg = sparse(Iidx,Jidx, Vidx, pp.n_elem, pp.n_elem );

% find elems which are connected to elems ee
function elems= find_adjoin(ee, ELEM)
   nn= ELEM(:,ee);
   [d,e]= size(ELEM);
   ss= zeros(1,e);
   for i=1:d
     ss= ss+ any(ELEM==nn(i));
   end
   elems= find(ss==d-1);
