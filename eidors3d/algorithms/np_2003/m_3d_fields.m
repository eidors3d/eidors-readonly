function [v_f] = m_3d_fields(vtx,el_no,m_ind,E,tol,gnd_ind,v_f);
%function [v_f] = m_3d_fields(vtx,el_no,m_ind,E,tol,gnd_ind,v_f);
%
%This function calculates the measurement fields using preconditioned conjugate gradients.
%
%
%
%vtx     = The vertices
%el_no   = The total number of electrodes in the system
%m_ind   = The measurements matrix (indices of electrode pairs)
%E       = The full rank system matrix
%tol     = The tolerance in the forward solution 
%gnd_ind = The ground index
%v_f     = The measurements fields

[vr,vc] = size(vtx);

Is_supl = zeros(vr,size(m_ind,1)); 
%no of electrodes x no of measurements (now currents)!

MC = [];

for i=1:size(m_ind,1)
   
   m_n = zeros(el_no,1);
   
   m_n(m_ind(i,1)) = 1;
   m_n(m_ind(i,2)) = -1;
   
   MC = [MC,m_n];
   
end

I = [Is_supl;MC];
I(gnd_ind,:) = 0;

if nargin < 7
v_f = zeros(size(E,1),size(I,2));
end

if exist('OCTAVE_VERSION')
   v_f= E\I;
   return;
end



if isreal(E)==1

%Preconditioner
M = cholinc(E,tol/10);

    for y=1:size(MC,2)
    %Set this line to suit your approximation needs. ***************
    %for more details use help pcg on Matlab's command window.
    [v_f(:,y),flag,relres,iter,resvec] = pcg(E,I(:,y),tol*norm(I(:,y)),1000,M',M,v_f(:,y)); 
    end
end %is real


if isreal(E)==0

[L,U] = luinc(E,tol/10);

  for y=1:size(MC,2)

  [v_f(:,y),flag,relres,iter,resvec] = bicgstab(E,I(:,y),tol*norm(I(:,y)),1000,L,U);

  end 
end %is complex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
