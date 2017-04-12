space_list = linspace(0.2, 0.03, 6.0);
angl = [330,200,80,100];
for ii = 1:length(space_list);
   space = space_list(ii);
   ypts = space*[1,-1];
   ytxt = sprintf('sphere(0.0,%f,0;0.01) or ',ypts); ytxt = ytxt(1:end-3);
   extra={'els',['solid els = ', ytxt, ' -maxh=0.01;']};
   fmdl= ng_mk_cyl_models([0,1,.15],angl(:),[0.1,0,0.02],extra);
   xnode = fmdl.nodes(:,1); ynode = fmdl.nodes(:,2);
   for i = 1:length(ypts);
      fmdl.electrode(end+1).nodes = find( (xnode.^2 + (ynode-ypts(i)).^2 <= .011^2));
      fmdl.electrode(end  ).z_contact = fmdl.electrode(1).z_contact;
   end
   img(ii) = mk_image(fmdl,1);
   clf; show_fem(img(ii),[0 1]); axis off;
   print_convert(sprintf('model_reduction01%s.png','a'+ii-1));
end
