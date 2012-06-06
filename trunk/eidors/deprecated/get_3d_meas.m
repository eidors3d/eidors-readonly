function [voltageH,voltageV,indH,indV,df] = get_3d_meas(elec,vtx,V,Ib,no_pl);
% GET_3D_MEAS: extracts multiplane voltage measurements from a calculated
% 3D nodal potential distribution V inside a tank with (no_pl) electrode
% planes. Each plane holds the same number of electrodes. Only the
% non-current carring electrodes at the time are involved in the
% measurements.
%
% [voltageH,voltageV,indH,indV,df]=get_3d_meas(elec,vtx,V,Ib,no_pl);
%
%elec      = The electrodes matrix
%vtx       = The vertices
%V         = The calculated forward solution
%Ib        = The current patterns without the zeroes patch
%no_pl     = The number of planes
%
%voltage_H = Horrisontal (local plane) measurements
%indH      = The two column matrix indicating the indices of the electrodes
%            involved in voltage_H, e.g. indH = [2 3; 3 4; 4 5;...] implies
%            voltage_H(1) = voltage(2) - volatge(3), etc
%df        = Array indexing the (numbero of) measurements to their corresponding
%				 current patterns.
%voltage_V = Vertical interplanar measurements
%indV      = ...


warning('EIDORS:deprecated','GET_3D_MEAS is deprecated as of 06-Jun-2012. ');

if size(V,2)~= size(Ib,2)
   error('Unmatched pattens')
end

[el_no,q] = size(elec);

el_pp = el_no/no_pl;

a=1:el_no;

X = reshape(a,el_pp,no_pl)';


Vm = V(size(vtx,1)+1:size(V,1),:); %Lower chunk of forward solution (complete electrode model)

voltageH = [];
indH = [];

df = [];

for w=1:size(Vm,2) %For each column of Vm

   cn = 0; %RESET the count of measurements per injection

   this_inj = Vm(:,w); %(no_of_electrodes x 1) vector

   for vv = 1:el_pp:el_no %i.e. 1 17 33 49 for 4 planes of 16 electrodes

      for t=vv:vv+(el_pp-1)-1 %t=1:15

         if Ib(t,w) == 0  && Ib(t+1,w) == 0   %Electrode not in the drive pair

            voltageH = [voltageH; (this_inj(t)-this_inj(t+1))];
            indH = [indH;[t , t+1]];
            cn = cn+1;
         end

         if t == vv+(el_pp-1)-1 && Ib(vv,w) == 0 && Ib(t+1,w) == 0

            voltageH = [voltageH; (this_inj(t+1))-this_inj(vv)]; %or is it vv=1;
            indH = [indH;[t+1, vv]];
            cn = cn+1;
         end

      end %for t -Measurements of the one plane

   end %for vv -Measurements for all electrode planes

   df = [df;cn];

   voltageV = [];
   indV = [];

   Y = reshape(X,el_no,1);


   cn = 0;
   wc = w;


   this_inj = Vm(:,wc); %(no_of_electrodes x 1) vector

   for ee = 1:no_pl:el_no

      this_chunk = Y(ee:ee+no_pl-1);

      for jj=1:length(this_chunk)-1

         if Ib(this_chunk(jj),wc) == 0 && Ib(this_chunk(jj+1),wc) == 0 %Electrodes not involved in currents

            voltageV = [voltageV; ((this_inj(this_chunk(jj)))- this_inj(this_chunk(jj+1)))];
            indV = [indV;[this_chunk(jj),this_chunk(jj+1)]];
            cn = cn+1;
         end

      end

   end
   df = [df;cn];
end


%voltage = [voltageH;voltageV];

%ind = [indH;indV];

% Separate df (Horrizontal / Vertical electrode combinations per current pattern as)
% dfh = df(1:2:end);
% dfv = df(2:2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

