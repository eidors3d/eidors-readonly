function [Vfields] = m_2d_fields(vtx,elec,E,pp,gnd_ind)
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

[er,ec] = size(elec);
[vr,vc] = size(vtx);
%First we need to set up the fields by setting up the relevant currents
if protocol == '{op}'
   d = er./2;
   df = [1:d];
   ds = [d+1:2*d];
   p = [];
   
for i=1:d
     
   for k=1:er
      
      if k ~= df(i) & k~= ds(i)
         p = [p;k];
      end
   end
end

p_p = reshape(p,length(p)./d,d); %Each column for each current exc.

ppp = [tril(p_p);triu(p_p,1)];

pc = ppp(:);

ps = nonzeros(pc);

mc = [];

for j=1:length(ps)-1
   
   this_mc = zeros(er,1);
   
   dop = ps(j); don = ps(j+1);
   
   if don - dop == 1 | dop - don == er-1
       
   this_mc(dop) = 1; this_mc(don) = -1;
     
   mc = [mc,this_mc];
   end
   
end

if model =='{ful}'
Is_supl = zeros(vr,(er*(er-4)/2));
I = [Is_supl;mc];
I(gnd_ind,:) = 0;

I = I(pp,:);
[Vfields] = forward_solver(vtx,E,I);
rr(pp) = 1:max(size(pp));
Vfields = Vfields(rr,:);
end

if model == '{gap}'
I = [];
  for h=1:size(mc,2)
      Ip = zeros(length(vtx),1);
      t_c = mc(:,h); %mc's hth column
      for tt=1:length(t_c)
         if t_c(tt) ~= 0
	        for g=1:size(elec,2)
               Ip(elec(tt,g)) = t_c(tt)./size(elec,2);
            end
         end
      end
   I = [I,Ip];
  end
  I(gnd_ind,:) = 0;
I = I(pp,:);
[Vfields] = forward_solver(vtx,E,I);
rr(pp) = 1:max(size(pp));
Vfields = Vfields(rr,:);
end

end % 'op'

%----------------------------------------------------------------------
if protocol == '{ad}'
   d = er;
   df = [1:d];
   ds = [2:d,1];
   p = [];
   
for i=1:d
     
   for k=1:er
      
      if k ~= df(i) & k~= ds(i)
         p = [p;k];
      end
   end
end

p_p = reshape(p,length(p)./d,d); %Each column for each current exc.

ppp = [tril(p_p);triu(p_p,1)];

pc = ppp(:);

ps = nonzeros(pc);

mc = [];

for j=1:length(ps)-1
   
   this_mc = zeros(er,1);
   
   dop = ps(j); don = ps(j+1);
   
   if don - dop == 1 | dop - don == er-1
       
   this_mc(dop) = 1; this_mc(don) = -1;
     
   mc = [mc,this_mc];
   end
   
end

if model =='{ful}'
Is_supl = zeros(vr,(er*(er-3)));
I = [Is_supl;mc];
I(gnd_ind,:) = 0;
I = I(pp,:);
[Vfields] = forward_solver(vtx,E,I);
rr(pp) = 1:max(size(pp));
Vfields = Vfields(rr,:);
end

if model == '{gap}'
   
I = [];
   for h=1:size(mc,2)
      Ip = zeros(length(vtx),1);
      t_c = mc(:,h); %mc's hth column
      for tt=1:length(t_c)
         if t_c(tt) ~= 0
            for g=1:size(elec,2)
               Ip(elec(tt,g)) = t_c(tt)./size(elec,2);
            end
         end
      end
   I = [I,Ip];
  end
 
I(gnd_ind,:) = 0;
I = I(pp,:);
[Vfields] = forward_solver(vtx,E,I);
rr(pp) = 1:max(size(pp));
Vfields = Vfields(rr,:);
end

end %'ad'

%------------------------------------------------------------------------
if protocol == '{co}'
   d = er./2;
	df1 =[1:d];
	df2 =[d+1:er];
	df3 =[2:d+1];
	df4 =[d+2:er,1];

   p = [];
   
for i=1:d
     
   for k=1:er
      
      if k ~= df1(i) & k ~= df2(i) & k ~= df3(i) & k ~= df4(i)
         p = [p;k];
      end
   end
end

p_p = reshape(p,length(p)./d,d); %Each column for each current exc.

pp1 = p_p(:,1:end-1);

ppp = [tril(pp1);triu(pp1,1)];

pc = ppp(:);

ps1 = nonzeros(pc);

pp2 = p_p(:,end); 

ps2 = [pp2((length(pp2)/2)+1:end);pp2(1:length(pp2)/2)];

ps = [ps1;ps2];

mc = [];

for j=1:length(ps)-1
   
   this_mc = zeros(er,1);
   
   dop = ps(j); don = ps(j+1);
   
   if don - dop == 1 | dop - don == er-1
       
   this_mc(dop) = 1; this_mc(don) = -1;
     
   mc = [mc,this_mc];
   end
   
end

if model =='{ful}'
Is_supl = zeros(vr,(d*(er-6)));
I = [Is_supl;mc];
I(gnd_ind,:) = 0;
I = I(pp,:);
[Vfields] = forward_solver(vtx,E,I);
rr(pp) = 1:max(size(pp));
Vfields = Vfields(rr,:);
end

if model == '{gap}'
   
I = [];
   for h=1:size(mc,2)
      Ip = zeros(length(vtx),1);
      t_c = mc(:,h); %mc's hth column
      for tt=1:length(t_c)
         if t_c(tt) ~= 0
            for g=1:size(elec,2)
               Ip(elec(tt,g)) = t_c(tt)./size(elec,2);
            end
         end
      end
   I = [I,Ip];
  end
 
I(gnd_ind,:) = 0;
I = I(pp,:);
[Vfields] = forward_solver(vtx,E,I);
rr(pp) = 1:max(size(pp));
Vfields = Vfields(rr,:);

end

end % '{co}'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KURSK Ver 2 							  								%
% Developed by: Nick Polydorides 								%
% First year Ph.D. achievement									%
% Copyright (c) N. Polydorides, September 2000			   %	
% Required: MATLAB 5.3 or update									%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
