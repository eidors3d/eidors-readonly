clear fmdl;
fmdl.name= 'Thorax FEM';
fmdl.type= 'fwd_model';
fmdl.elems= tri;
fmdl.nodes= vtx;
fmdl.boundary= double( find_boundary(tri) );

k=1;
for i= 1:size(fmdl.boundary,1)
   fmbi= fmdl.boundary(i,:);
   if any(fmbi(1) == enodes) &&  ...
      any(fmbi(2) == enodes) &&  ...
      any(fmbi(3) == enodes); 
       xy_elec= mean( vtx(fmbi,1:2) );
       angl= atan2(xy_elec(1),xy_elec(2));
       ak(k)= (angl);
%      fmdl.electrode(k).nodes = fmbi;
       k=k+1;
    end
end
[n,x]= hist(ak,100);
up = find(diff(n==0) ==1)+1;
dn = find(diff(n==0) ==-1);
sep = [-10,mean([up;dn],1),110];
plot(n,'-');hold on;plot(sep,1,'k*');hold off
sep = ((sep-1)/99 - .5)*2*pi;

% female
sep=[ ...
 -2.8560, -2.4415, -2.1036, -1.7505, -1.4587, -1.1209, -0.7370, -0.3378, ...
  0.2150,  0.6756,  1.0749,  1.4127,  1.7658,  2.1497,  2.5489,  2.9942];

%male
sep=[ ...
 -2.8560, -2.3340, -1.8887, -1.4587, -1.1516, -0.8906, -0.5681, -0.1996, ...
  0.1843,  0.5681,  0.8599,  1.1823,  1.5355,  1.9040,  2.3186,  2.7946]


eangl= atan2(vtx(enodes,1),vtx(enodes,2));
for i=1:16
   if i==1; 
      this_elec= find(eangl <  sep(1) | eangl >= sep(16)); 
   else
      this_elec= find(eangl >= sep(i-1) & eangl < sep(i)); 
   end
   fmdl.electrode(i).nodes = enodes(this_elec);
   fmdl.electrode(i).z_contact= 0.1;
end
fmdl.stimulation = [];
fmdl.gnd_node= 1;
show_fem(fmdl);
