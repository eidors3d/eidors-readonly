function [JTb] = adjoint_spin(vtx,simp,elec,x,gnd_ind,zc,I,no_pl,Vmes);
%function [JTb] = adjoint_spin(vtx,simp,elec,x,gnd_ind,zc,I,no_pl,Vmes);
%
%The function calculates the product J'*b, i.e. Jacobian transpose times 
%a (measurements Vmes) vector b, using the adjoint sources formulation.
% 
%
%
%JTb     = J'*Vmes
%vtx     = The vertices matrix
%simp    = The simplices matrix
%elec    = The electrodes matrix
%x       = The solution estimate based on which J is supposed to be caclulated
%gnd_ind = The grounf index
%zc      = The electrode contact impedance
%I       = The current patterns
%no_pl   = Number of electrode planes
%Vmes    = The measurements vector


spin = 1;%size(I,2);
[vr] = size(vtx,1);
ptr = 0;
    
     %Update model based on last x 
     [E,D,Ela,pp] = fem_master_full(vtx,simp,x,gnd_ind,elec,zc,'{n}');
     EB = kron(eye(spin),inv(E));    
     
          
     for j=1:spin:size(I,2); %For each group of current patterns
     
     Ib = I(:,j:j+spin-1);
     IB = reshape(Ib,size(Ib,1)*size(Ib,2),1);
     
     VB = EB*IB;
         
     
     %The current fields for I(: ,spin) only
     VB1 = reshape(VB,size(vtx,1)+size(elec,1),spin);
     VB2 = VB1(1:size(vtx,1),:); 
     %These are the standard current fields used in the calculation of the Jacobian
     Vjcf = reshape(VB2,size(vtx,1)*spin,1);
     %Electrode potentials removed
     
     volts = []; %Some simulated data based on x 
     ind = [];
     dfh = [];
     
     [voltageH,voltageV,indH,indV,df] = get_3d_meas(elec,vtx,VB1,... 
                                                       I(size(vtx,1)+1:end,j:j+spin-1),no_pl);
     volts = [volts;voltageH];
     ind = [ind;indH];
     dfh = [dfh;df(1:2:end)];
       
  
          
         %The measurement fields for the measurements during I(:,spin) are active
         VmI = Vmes(ptr+1:ptr+sum(dfh)); %Part of the measurements
         ptr  = ptr+sum(dfh);
     
                
         %Set up the measurement patterns for all I(:,spin)
         Imb = [];
         MC = [];
         kap = 1;

             for ij=1:size(ind,1)
                 m_n = zeros(size(elec,1),1);
                 m_n(ind(ij,1)) = 1;
                 m_n(ind(ij,2)) = -1;
                 MC = [MC,m_n];
                 if ij == sum(dfh(1:kap))
                     if ij ~= size(ind,1)
                     kap = kap+1;
                     end
                     Imb = sparse(blkdiag(Imb,MC(:,1+ij-dfh(kap):end)));
                 end
            end
             
                 
         %These are measurement residuals of the j'th current pattern
         b = VmI;
              
         Im = Imb*b; 
         
         %Reshape it in columns. 
         Im_col = reshape(Im,size(elec,1),spin);
         mul = 1;
         
                    I_s = [];
                    for t=1:spin
                        I_s = [I_s; zeros(size(vtx,1),1); Im_col(:,t)*mul];
                    end
     
                    Vmf = EB*I_s;
                    Vmf1 = reshape(Vmf,size(vtx,1)+size(elec,1),spin);
                    Vmf2 = Vmf1(1:size(vtx,1),:);
                    Vjmf = reshape(Vmf2,size(vtx,1)*spin,1); %
                    
                    DM = kron(eye(spin),D); % Gradient operator
                    Vjmf_n = -DM * Vjmf;
                    Vjcf_n = -DM * Vjcf;
                    JTb_p = Vjmf_n .* Vjcf_n ;  
                    J_m = reshape(JTb_p,size(simp,1)*3,spin);
                    if spin > 1
                    JTb_x3 = sum(J_m.');
                else
                    JTb_x3 = J_m.';
                end
                    JTb_u = JTb_x3(1:3:end) + JTb_x3(2:3:end) + JTb_x3(3:3:end);
        
                    JTb = -JTb_u.' .* diag(Ela(1:3:end,1:3:end));
        
   end
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   