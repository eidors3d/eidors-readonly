function [img,map]= GREIT_NOSER_ndiff( ref_meas, reconst_meas )
% Reconstruct GREIT images using NOSER algorithm
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   [RM,map] = calc_RM(1);

   % Expand ref_meas to the full size of reconst_meas
   num_meas = size(reconst_meas,2);
   ref_meas = ref_meas * ones(1,num_meas);
   dv = ( reconst_meas - ref_meas ) ./ ref_meas; % CHANGE IS HERE:

   % reconst image
   ds = RM*dv;

   img= reshape(ds, 32,32,num_meas);

function [RM,map] = calc_RM(noser_p, diff_p)
   [J,map,vbkgnd] = GREIT_Jacobian_cyl;
   J = J ./ (vbkgnd*ones(1,size(J,2))); % Normalized Jacobian
   RM = zeros(size(J'));

   L = mk_diff_prior(map);
 
   % Remove space outside FEM model
   J= J(:,map);
   % inefficient code - but for clarity
   diagJtJ = diag(J'*J).^(.25);
   D= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

%  R= D'*D; %NOSER
   R= D'*L'*L*D; %Diff
   hp = 0.1;
   RM(map,:)= (J'*J + hp^2*R)\J';

function L= mk_diff_prior(map)
   sz= prod(size(map));
   cube0 = ones(size(map));

   cube= cube0; cube(end,:) = 0; f1= find(cube(:)); lf1 = length(f1);
   cube= cube0; cube(1,:) =   0; f2= find(cube(:)); lf2 = length(f2);
   oo= cube(f1) + cube(f2);
   Lx= sparse(1:lf1,f1, oo,lf1,sz) + sparse(1:lf2,f2,-oo,lf2,sz);

   cube= cube0; cube(:,end) = 0; f1= find(cube(:)); lf1 = length(f1);
   cube= cube0; cube(:,1)   = 0; f2= find(cube(:)); lf2 = length(f2);
   oo= cube(f1) + cube(f2);
   Ly= sparse(1:lf1,f1, oo,lf1,sz) + sparse(1:lf2,f2,-oo,lf2,sz);

   L= [Lx;Ly];
   L= L(:,map);
