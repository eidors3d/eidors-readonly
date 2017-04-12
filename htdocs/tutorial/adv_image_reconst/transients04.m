%Create 2D model
ball1 = 'solid ball = cylinder(0.1,0.4,0;0.1,0.4,1;0.3) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;';
box1  = 'solid box = orthobrick(-0.6,-0.6,0; 0.1,-0.1,0.05) -maxh=0.1;';
fmdl= ng_mk_cyl_models(0,[8],[0.2,0,0.05],{'ball','box',[ball1,box1]}); 
fmdl.stimulation = stim_meas_list([1,5,2,3],8,.01,1);

%Conductivity and permittivity parameters
img = mk_image( fmdl, 1);
img.elem_data( fmdl.mat_idx{2} ) = .01 + .01i;
img.elem_data( fmdl.mat_idx{3} ) = .01 + .01i;
img.calc_colours.component = 'real';

subplot(221); show_fem(img);
print_convert transients01a.png


return
tlim = 1e+0; tpts = 300;
fmax = 1/tlim*(tpts/2);
 tax = linspace(0,tlim,tpts); % 1 sec
puls= tax> 0.4*tlim & tax<0.6*tlim;
fpul = fftshift(fft(puls));
fax = linspace(0,fmax,tpts+1); fax(end)=[];
fax = fftshift(fax);
fax(fax>0.5*fmax) = fax(fax>0.5*fmax) - fmax;
subplot(423); plot(tax, puls); ylim([-.1,1.1]);
subplot(424); plot(fax, real(fpul),'.');% ylim([-.1,1.1]);

red = real(img.elem_data);
ied = imag(img.elem_data);
img.fwd_solve.get_all_meas=1; %Get all measurements
clear vn; for i=fliplr(1:tpts)
   img.elem_data= red + ied*fax(i)*2j*pi;
   vv = fwd_solve(img);
   vn(:,i) = fpul(i)*vv.volt;
end

imgn = rmfield(img,'elem_data');
imgn.node_data = real( vn(:,floor(tpts*.5)+1  ));
subplot(425);  show_fem(imgn);
imgn.node_data = real( vn(:,floor(tpts*.5)+25 ));
subplot(426);  show_fem(imgn);

vt = ifft(ifftshift(vn,2),[],2);
if norm(real(vt(:))) < 1e-10*norm(imag(vt(:))); error('ifft'); end

subplot(414); 
imgn.node_data = real(vt);
imgn.show_slices.img_cols = 20;
rimg = calc_slices(imgn);

%show_slices( imgn);
 slice = rimg(:,32,:);
%slice = rimg(32,:,:);
image(calc_colours(reshape(slice,64,[]), imgn))
 xlim([.30,.70]*tpts);

