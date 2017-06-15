img.calc_colours.component = 'real'; % display real component of current

for i=1:2
   if i==1; cplx= .01;
   else   ; cplx= .01 - pi*2i *100; end
   img.elem_data( inclusion ) = cplx;

   subplot(2,2,i); show_current(img); axis off
end
print_convert transients02a.png

red =  real(img.elem_data); % in-phase
ied = (imag(img.elem_data)~=0)*.01; % out-of-phase
