% Test area of electrodes
% $Id: build_single_plane02.m,v 1.3 2007-09-24 14:29:56 aadler Exp $

% verify np_fwd_parameters identifies
for i= 1:electrodes_per_plane*number_of_planes
   elec_areas= [];
      fprintf('\nelec#%02d (area,nodes):',i);
   for fmdl= { ng_mdl_16x1_coarse, ng_mdl_16x1_fine, ng_mdl_16x1_vfine }
      pp= np_fwd_parameters( fmdl{1} );
      ee= reshape( pp.elec(i,:), 3, []);
      ee(:, all(ee==0,1)) = [];

   %  Xs = reshape(pp.vtx(ee,1),3,[]);
   %  Ys = reshape(pp.vtx(ee,2),3,[]);
   %  Zs = reshape(pp.vtx(ee,3),3,[]);
   %  h=patch(Xs,Ys,Zs, 'b');
   %  set(h, 'FaceLighting','none', 'CDataMapping', 'direct' );

      elec_area= 0; elec_nodes= size(ee,2);
      for j= 1:elec_nodes
         elec_area = elec_area + triarea3d( pp.vtx( ee(:,j),:) );
      end
      fprintf('[ %6.5f (%3d) ]',elec_area, elec_nodes);

   end
end

