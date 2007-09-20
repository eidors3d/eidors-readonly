% Test area of electrodes
% $Id: build_single_plane02.m,v 1.1 2007-09-20 20:30:33 aadler Exp $

% verify np_fwd_parameters identifies
pp= np_fwd_parameters( fmdl);
for i= 1:electrodes_per_plane*number_of_planes
   ee= reshape( pp.elec(i,:), 3, []);
   ee(:, all(ee==0,1)) = [];

   Xs = reshape(fmdl.nodes(ee,1),3,[]);
   Ys = reshape(fmdl.nodes(ee,2),3,[]);
   Zs = reshape(fmdl.nodes(ee,3),3,[]);
   h=patch(Xs,Ys,Zs, 'b');
   set(h, 'FaceLighting','none', 'CDataMapping', 'direct' );

   elec_area= 0;
   for j= 1: size(ee,2)
      elec_area = elec_area + triarea3d( fmdl.nodes( ee(:,j),:) );
   end

   fprintf('elec#%d: area=%f\n',i,elec_area);
end

show_fem( fmdl);
print -r100 -dpng build_single_plane01a.png;

