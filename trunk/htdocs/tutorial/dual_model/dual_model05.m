% Simulate data $Id$

% Reconstruct
for model= 1:3
   img(model)= inv_solve(imdl(model), vh, vi);
end

% Show image mapped to fine model
for model= 1:3
   subplot(2,3,model)
   show_fem(img(model));
end

print -r125 -dpng dual_model05b.png;

% Show image mapped to coarse model
for model= 1:3
   subplot(2,3,model)
   img(model).fwd_model = mdl_coarse.fwd_model;
   show_fem(img(model));
end

print -r125 -dpng dual_model05a.png;
