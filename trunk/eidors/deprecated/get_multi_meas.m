function [voltage,ind,df] = get_multi_meas(protocol,elec,V,I,vtx,no_pl);
%function [voltage,ind,df] = get_multi_meas(protocol,elec,V,I,vtx,no_pl);
%
%The function can be used in the occasions where plane current patterns
%(adjacent or polar) are adopted for systems with more planes, i.e. set by
%set_multi_currents function. Only non-current carrying electrodes are
%involved in the measurements.
%
%
%
%protocol= The selected protocol '{op}' or '{ad}'
%elec    = The electrodes
%no_pl   = The number of planes
%I       = The RHS vectors, i.e., the current patterns padded with zeroes
%          for the forward calculations
%V       = The forward solution
%voltage = The array of measureemnts according to the selected protocol
%ind     = The index of the measurements


warning('EIDORS:deprecated','GET_MULTI_MEAS is deprecated as of 06-Jun-2012. ');

if size(V,2)~= size(I,2)
   error('Unmatched pattens')
end


no_el = size(elec,1);
el_pp = no_el/no_pl;

I = I(size(vtx,1)+1:end,:);  %Lower chunk
Vm = V(size(vtx,1)+1:end,:); %Lower chunk


voltage = [];
ind = [];
df = [];

d = size(I,2); %Injections

for i=1:d %for each injection

   cn = 0;

   L = [];
   fst = 0;
   kk=0;

   for ej=1:el_pp:no_el

      kk=ej-1;

      for ei=1:el_pp-1
         if ei==1
            fst = ej;
         end
         L = [L;[kk+ei kk+(ei+1)]];
      end
      L = [L;[fst fst+16-1]];
   end

   LMP = []; %Exclude measurements engaging electrodes from different electrodes

   for j=1:size(L,1)

      if ceil(L(j,1)/el_pp) == ceil(L(j,2)/el_pp)

         LMP = [LMP;L(j,:)];
      end
   end

   el_in = find(I(:,i)~=0); % The current carring electrodes

   LM =[];

   for k=1:size(LMP,1)
      if ismember(LMP(k,:),el_in) == [0 0]
         LM = [LM;LMP(k,:)]; %Modified list current carring electrodes removed.
      end
   end

   lf = size(LM,1);
   dd = lf/no_pl;

   LN = [];

   k=0;

   for in=1:dd:lf %   1    13    25    37

      k = k+1;

      th_LM = LM(in:in+dd-1,:); %size(th_LM) = (12 2)

      if protocol == '{op}'
         for j=1:size(th_LM,1)
            if th_LM(j,1)== intersect(find(I(:,i)>0),in:k*el_pp)+1
               LN = [th_LM(j:end,:);th_LM(1:j-1,:)];
               break;
            end
         end
      elseif protocol =='{ad}'
         for j=1:size(th_LM,1)
            if th_LM(j,1)== intersect(find(I(:,i)<0),in:k*el_pp)+1
               LN = [th_LM(j:end,:);th_LM(1:j-1,:)];
               break;
            end
         end
      else
         error('Protocol needs to be {op} or {ad}');
      end

      for u=1:size(LN,1)
         voltage = [voltage; (Vm(LN(u,1),i)- Vm(LN(u,2),i))];
      end


      ind = [ind;LN];
      cn = cn+1;

   end %for in

   df = [df;cn];
end %for injections





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




