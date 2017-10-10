%MIKE
%0. Compare free space and Neumann greens
%1. Use true Dz (N(x,z) in phess_obj
%2. Find error in Taylor Hessian - DONE with test_saturation_hessian.. 
%showing (some) validation of this..
%3. Make objective Hessian with -  DONE testing in test_satuation_hessian by
%plotting the DF'DF hessian and H*DV hessian as function of perturbation
%which demonstrates that as pertrubation decreases hessian of obejctive is
%converging to GN only hessian. BUT needs further valdiation AND testing 
%with phess (see FRISS 1 below)


%FRISS
%0. Change ellement select in first order expansion i.e. optimise over
%interior region + make a mesh
%1. TaylorPhessian of obejctive comparison
%2. Inverse solve with BFGS ..
%3. Compare BFGS Hessian to Taylor Hessina

%NEWTODO
%1 convex vs non-convex problem
%2 long triangle - deviatoric part of M tensor ->> 0 more rapidly?
%3 fix sensitivity of p-tensor and compare, incl. rotated/shaped els
%4 reconstruction model which doesn't include els at bdry
%
% take home - less accurate but works for large scale?