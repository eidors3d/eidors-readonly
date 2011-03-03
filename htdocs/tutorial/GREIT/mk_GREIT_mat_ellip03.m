% $Id$

for k=1:4
   rimg = inv_solve(imdl(k),vh,vi); %Reconstruct

   subplot(1,4,k)
   show_fem(rimg); axis square;
end

print_convert mk_GREIT_matrix03.png '-density 180'
