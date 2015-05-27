function imdl_2x16 = get_base_model

if exist('imdl_2x16.mat','file');
   load imdl_2x16;
   return
else
   [fmdl_2x16, midx] = ng_mk_cyl_models([0.38,0.145,0.0055],[16,0.19,0.262],[0.001, 0, 0.005]);
   %% re-position electrodes
   % even electrodes are below odd ones
   idx = [2:2:32 1:2:32 ];
   elecs = fmdl_2x16.electrode;
   for i = 1:32;
      fmdl_2x16.electrode(idx(i)) = elecs(i);
   end
   %%
   imdl_2x16 = select_imdl(fmdl_2x16,{'Basic GN dif'});
   opt.imgsz = [32 32];
   opt.zvec = [-.001 .14:.02:.32 .381];
   opt.prune_model= false;
   imdl_2x16 = mk_voxel_volume(imdl_2x16,opt);
   % prune manually
   rmdl = imdl_2x16.rec_model;
   D = sqrt(sum(rmdl.nodes(:,1:2).^2,2)) <=.145;
   E = all(D(reshape(rmdl.elems',4*6,[])'),2);
   E2= logical(kron(E,true(6,1)));
   rmdl.elems = rmdl.elems(E2,:);
   rmdl.coarse2fine = rmdl.coarse2fine(E2,E);
   rmdl = rmfield(rmdl, 'boundary');
   rmdl.inside = E2;
   imdl_2x16.rec_model = rmdl;
   save imdl_2x16 imdl_2x16

end
