skip_list=[0 2]; % skip N electrodes, stim/meas pattern
ne=8; % number of electrodes
ns=5; % number of stream lines
clf;
for ii=1:length(skip_list)
   ii
   skip_elec=skip_list(ii)
   imdl = mk_common_model('e2C',ne);
   stim=mk_stim_patterns(ne,1,[1 2+skip_elec],[1 2+skip_elec],{},1);
   stim=stim_meas_list(stim);
   stim=sortrows(stim,[1 2 4 3]);
   imdl.fwd_model.stimulation = stim_meas_list(stim);
   img=mk_image(imdl,1);
   img.fwd_solve.get_all_meas=1; % data.volt = all FEM nodes, but not CEM
   data=fwd_solve(img);
   img=rmfield(img,'elem_data');
   for sel = [1 2 3];
      switch sel
      case 1
         stimi = 1; measj = 1; ij=1;
      case 2
         stimi = 1; measj = 2; ij=2;
      case 3
         stimi = 2; measj = ii; ij=(ne-3)+ii;
      otherwise
         error('oops');
      end
      stim(ij,:)
      nn = imdl.fwd_model.nodes;
      % stim streamlines
      imdl.fwd_model.stimulation = stim_meas_list(stim(ij,:),ne);
      imgt=mk_image(imdl,1);
      imgt.fwd_solve.get_all_meas=1; % data.volt = all FEM nodes, but not CEM
      datat=fwd_solve(imgt);
      PLANE= [inf,inf,0]; % show voltages on this slice
      imgt.fwd_model.mdl_slice_mapper.npx = 64;
      imgt.fwd_model.mdl_slice_mapper.npy = 64;
      imgt.fwd_model.mdl_slice_mapper.level = PLANE;
      q = show_current(imgt, datat.volt(:,1));
      clear sxy;
      for jj=[1 2]
         s=stim(ij,jj);
         ee = imdl.fwd_model.electrode(s).nodes;
         sxy(jj,:) = (max(nn(ee,:),[],1)+min(nn(ee,:),[],1))./2;
      end
      ssxy=mean(sxy); theta=atan2(ssxy(1),ssxy(2))+pi/2;
      R = [ cos(theta) -sin(theta); sin(theta) cos(theta) ];
      sx = linspace(-.5,.5,ns)';
      sxy = [sx sx*0] * R;
      sx = sxy(:,1); sy = sxy(:,2);
      % measj streamlines
      meas = stim(:,[3 4 1 2]);
      imdl.fwd_model.stimulation = stim_meas_list(meas(ij,:),ne);
      imgt=mk_image(imdl,1);
      imgt.fwd_solve.get_all_meas=1; % data.volt = all FEM nodes, but not CEM
      datat=fwd_solve(imgt);
      imgt.fwd_model.mdl_slice_mapper.npx = 64;
      imgt.fwd_model.mdl_slice_mapper.npy = 64;
      imgt.fwd_model.mdl_slice_mapper.level = PLANE;
      qq = show_current(imgt, datat.volt(:,1));
      clear mxy;
      for jj=[1 2]
         m=stim(ij,jj+2);
         ee = imdl.fwd_model.electrode(m).nodes;
         mxy(jj,:) = (max(nn(ee,:),[],1)+min(nn(ee,:),[],1))./2;
      end
      mmxy=mean(mxy); theta=atan2(mmxy(1),mmxy(2))+pi/2;
      R = [ cos(theta) -sin(theta); sin(theta) cos(theta) ];
      mx = linspace(-.5,.5,ns)';
      mxy = [mx mx*0] * R;
      mx = mxy(:,1); my = mxy(:,2);
      % find measj electrodes
      imdl.fwd_model.stimulation = stim_meas_list(stim,ne);
      mel=find(imdl.fwd_model.stimulation(stimi).meas_pattern(measj,:)~=0);
      mxy = [];
      for m=mel(:)'
         nn = imdl.fwd_model.nodes;
         ee = imdl.fwd_model.electrode(m).nodes;
         mxy(end+1,:) = (max(nn(ee,:),[],1)+min(nn(ee,:),[],1))./2;
      end
      img.node_data = data.volt(:,stimi);
      %subplot(2,3,(ii-1)*3+sel);
      clf; subplot(4,4,1);
      h=show_fem(img,[0 1]);
      set(h,'EdgeColor','none'); axis off;
      hold on; r = 0.12;
      rectangle('Position',[-1 -1 2 2],'Curvature',[1,1]);
      %rectangle('Position',[mxy(1,:)-r r*2 r*2],'Curvature',[1,1]);
      %rectangle('Position',[mxy(2,:)-r r*2 r*2],'Curvature',[1,1]);
      hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',2, 'color','b');
      hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',2, 'color','b');
      hh=streamline(qq.xp,qq.yp, qq.xc, qq.yc,mx,my); set(hh,'Linewidth',2, 'color','r');
      hh=streamline(qq.xp,qq.yp,-qq.xc,-qq.yc,mx,my); set(hh,'Linewidth',2, 'color','r');
      if 0
         hh=plot(sx,sy,'bo');
         hh=plot(mx,my,'ro');
         plot(ssxy(1),ssxy(2),'ko');
         plot(mmxy(1),mmxy(2),'ko');
      end

      hold off;
      print('-dpng',sprintf('skip%d%d.png',ii,sel));
   end
   stim_meas_list(imdl.fwd_model.stimulation)
end
