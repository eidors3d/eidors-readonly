imb=  mk_common_model('a2c2',16);
img= mk_image(imb,1); %show_fem(img,[0,0,1]);
vv = fwd_solve(img); vv = vv.meas;
p1 = [29,28]; p2 = [19,20];

sp = logspace(-2.5,2.5,40);
for i = length(sp):-1:1
   for j = length(sp):-1:1
      img.elem_data(p1) = sp(i);
      img.elem_data(p2) = sp(j);
      v2 = fwd_solve(img); 
      M(i,j) = norm(vv-v2.meas);
   end
end
img.elem_data(p1) = 10^(+.1);
img.elem_data(p2) = 10^(-.1);
subplot(121);show_fem(img,[0,0,1]);
clim = [min(M(:)),max(M(:))]*(eye(2)-[1,0;1,0]/20);
subplot(122);imagesc(log10(sp),log10(sp),M,clim); axis image

print_convert minima_region_01a.jpg
