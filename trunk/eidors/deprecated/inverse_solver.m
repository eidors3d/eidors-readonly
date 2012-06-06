function [solf,solp] = inverse_solver(I,voltage,tol,mat_ref,vtx,simp,elec,no_pl,zc,perm_sym,gnd_ind,tfac,Reg,it);
%function [solf,solp] = inverse_solver(I,voltage,tol,mat_ref,vtx,simp,elec,no_pl,zc,perm_sym,gnd_ind,tfac,Reg,it);
%
%Calculates a Newton non-linear inverse solution by iteration.
%
%
%
%solf    = The non-linear inverse solution
%mat_ref = Initial guess on the solution
%I       = The current patterns
%elec    = The electrode faces
%zc      = The contact impedance vector
%voltage = The measurements
%tfac    = The regularisation parameter
%Reg     = The regularisation matrix
%it      = Number of iterations
%vtx     = The vertices matrix
%simp    = The simplices matrix
%gnd_ind = The ground index (node/electrode)
%no_pl   = The number of planes


warning('EIDORS:deprecated','INVERSE_SOLVER is deprecated as of 06-Jun-2012. ');

tol = 1e-4; %Inverse calculations error tolerance. Change accodringly

sol_upd = mat_ref; %Initial estimate - homogeneous background

solp = zeros(size(simp,1),1); %Total change (plotting purposes)

Ib = I(size(vtx,1)+1:end,:);

el_no = size(elec,1);

if it==1
   disp('A linear solution will be calculated')
end


for i=1:it

   [E,D,Ela,pp] = fem_master_full(vtx,simp,sol_upd,gnd_ind,elec,zc,perm_sym);

   if i==1
      %sprintf('Current fields for iteration %d',i)
      [V] = forward_solver(E,I,tol,pp);
      [viH,viV,indH,indV,df] = get_3d_meas(elec,vtx,V,Ib,no_pl);
      dfv = df(1:2:end);
      vi = viH;
      %sprintf('Measurement fields for iteration %d',i)
      [v_f] = m_3d_fields(vtx,el_no,indH,E,tol,gnd_ind);
   else
      %sprintf('Current fields for iteration %d',i)
      [V] = forward_solver(E,I,tol,pp,V);
      [viH,viV,indH,indV,df] = get_3d_meas(elec,vtx,V,Ib,no_pl);
      dfv = df(1:2:end);
      vi = viH;
      %sprintf('Measurement fields for iteration %d',i)
      [v_f] = m_3d_fields(vtx,el_no,indH,E,tol,gnd_ind,v_f);
   end

   [J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,sol_upd,zc,v_f,dfv,tol,perm_sym);

   sol = (J.'*J + tfac*Reg.'*Reg)\ (J.' * (voltage - vi));
   %sol = pinv([J;sgrt(tfac)*Reg],tol) * [sqrt(tfac)*Reg*sol_upd; (voltage - vi)];
   sol_upd = sol_upd + sol;
   solp = solp + sol;

   h1 = figure;
   set(h1,'NumberTitle','off');
   set(h1,'Name','Reconstructed conductivity distribution');
   subplot(2,3,1); [fc] = slicer_plot_n(2.63,sol_upd,vtx,simp); view(2); grid; colorbar; axis off; title('z=2.63');
   subplot(2,3,2); [fc] = slicer_plot_n(2.10,sol_upd,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10');
   subplot(2,3,3); [fc] = slicer_plot_n(1.72,sol_upd,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72');
   subplot(2,3,4); [fc] = slicer_plot_n(1.10,sol_upd,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10');
   subplot(2,3,5); [fc] = slicer_plot_n(0.83,sol_upd,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
   subplot(2,3,6); [fc] = slicer_plot_n(0.10,sol_upd,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');
   drawnow;


   sprintf('Error norm at iteration %d is %f',i,norm(voltage - vi))

end %for it

solf = sol_upd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 6.1 R12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
