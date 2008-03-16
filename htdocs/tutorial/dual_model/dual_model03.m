% Simulate data $Id: dual_model03.m,v 1.2 2008-03-16 11:06:27 aadler Exp $

for model= 1:3
   img= inv_solve(imdl(model), vh, vi);
   subplot(2,3,model)
   show_fem(img);
end

print -r125 -dpng dual_model03a.png;
