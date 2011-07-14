subplot(221);
for loop1= 1:5;
       
  for loop2 = 1:2;
     if loop2 ==1; fn = sprintf('%d-control.RAW',loop1);
     else          fn = sprintf('%d-injury.RAW',loop1);
     end

     prname= sprintf('if_lung_injury02-%d%d.png', loop1,loop2);

     dd= eidors_readdata(fn);

     img=inv_solve(imdl,mean(dd,2),dd);
     [jnk,fmin] = min(mean(img.elem_data,1)); % find end-inspi
     [jnk,fmax] = max(mean(img.elem_data,1)); % find end-expi

     img=inv_solve(imdl,dd(:,fmax), dd(:,fmin));
     show_slices(img);

     print_convert(prname,'-density 50');
   end
end
