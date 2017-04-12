tlim =  0.05; tpts = 100;
fmax = 1/tlim*(tpts/2);
tax = linspace(0,tlim,tpts); 
puls= tax> 0.4*tlim & tax<0.6*tlim;
fpul = fftshift(fft(puls));
fax = linspace(0,fmax,tpts+1); fax(end)=[];
fax = fftshift(fax);
fax(fax>0.5*fmax) = fax(fax>0.5*fmax) - fmax;
subplot(221); plot(tax, puls); ylim([-.1,1.1]); box off;
subplot(222); plot(fax, real(fpul),'.');        box off;

clear vn; for i=fliplr(1:tpts)
   img.elem_data= red + ied*fax(i)*2j*pi;
   vv = fwd_solve(img);
   vn(:,i) = fpul(i)*vv.volt;
end
print_convert transients03a.png
