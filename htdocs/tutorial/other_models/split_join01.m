imdl = mk_common_model('a2c0',8); fmdl= imdl.fwd_model;
subplot(131); show_fem(fmdl ,[0,1,1]); xlim([-1,1]);

fmdl1= crop_model(fmdl,inline('x<-0.4','x','y','z'));
fmdl.nodes(:,1) = fmdl.nodes(:,1) + 0.25;
subplot(132); show_fem(fmdl1,[0,1,1]); xlim([-1,1]);

% reverse
fmdl2 = fmdl1; fmdl2.nodes(:,1) = -fmdl2.nodes(:,1);
fmdl2 = crop_model(fmdl2,inline('x+y<-1.28','x','y','z'));
fmdl2 = crop_model(fmdl2,inline('y<-0.95','x','y','z'));
fmdl2 = crop_model(fmdl2,inline('x<-0.95','x','y','z'));
subplot(133); show_fem(fmdl2,[0,1,1]); xlim([-1,1]);

print_convert('split_join01a.png',struct('pagesize',[12,5]));
