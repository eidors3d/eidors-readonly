% Inverse Transform
vt = ifft(ifftshift(vn,2),[],2);
if norm(real(vt(:))) < 1e-10*norm(imag(vt(:))); error('ifft'); end

subplot(211); 

imgn = rmfield(img,'elem_data');
imgn.node_data = real(vt);
imgn.show_slices.img_cols = 10;
imgn.get_img_data.frame_select = round(0.3*tpts):round(0.79*tpts);
show_slices( imgn);

print_convert transients04a.png
