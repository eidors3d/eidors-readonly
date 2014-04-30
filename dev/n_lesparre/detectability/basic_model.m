function img=basic_model

   [fmdl,mi1,mi2] = test_model(-1,-2,0.4,0.4);
   skip = [0,9];
   fmdl.stimulation = mk_stim_patterns(35,1,skip, skip, {},1);
   contrast = 1.1; N_sims = 1e3; SNR_noise = 1e5; SNR_ampl = 1e5; 

   img1= mk_image(fmdl,1); img1.elem_data(mi1) = contrast;
   v1 = do_simulations(img1, N_sims, SNR_noise, SNR_ampl);
   img2= mk_image(fmdl,1); img2.elem_data(mi2) = contrast;
   v2 = do_simulations(img2, N_sims, SNR_noise, SNR_ampl);
   subplot(221); show_fem(img1); view(0,0);
   subplot(222); show_fem(img2); view(0,0);

   u1 = mean(v1,2);
   u2 = mean(v2,2);
   C1 = spdiag( var( v1, 1, 2) );
   C2 = spdiag( var( v2, 1, 2) );
   w  = (C1 + C2)\(u1 - u2);
   w  = w/norm(w);
   w  = w/1000;
   subplot(223)
   [jnk,idx] = sort(abs(w));
   f = idx(end); g = idx(end-1);
   plot(v1(f,:),v1(g,:),'bo', v2(f,:),v2(g,:),'ro', v1(f) + [0;w(f)],v1(g) + [0;w(g)],'g-')

   d1 = w'*v1;
   d2 = w'*v2;

   tstar = (mean(d1) - mean(d2))/sqrt( (var(d1) + var(d2))/2 );

   subplot(224); cla;
   hist(d1); hold on; hist(d2); hold off; 
h = findobj(gca, 'Type','patch');
set(h(1), 'FaceColor','r', 'EdgeColor','k');
set(h(2), 'FaceColor','b', 'EdgeColor','k');
   title(sprintf('tstar = %5.3f',tstar));

   img = img1;


function vv = do_simulations(img, N_sims, SNR_noise, SNR_ampl)
   vi = fwd_solve(img);
   J  = calc_jacobian(img);
   vv = vi.meas*ones(1,N_sims) + ...
        randn(size(vi.meas,1),N_sims)/SNR_noise + ...
        J * randn(size(J,2),N_sims) / SNR_ampl;

   return
   imgs = img;
   for i= N_sims:-1:1
      imgs.elem_data = img.elem_data + randn(size(imgs.elem_data))/SNR_ampl;
      vi = fwd_solve(imgs);
      vv(:,i) = vi.meas + randn(size(vi.meas))/SNR_noise;
   end

function [fmdl,mi1,mi2] = test_model(xctr,zctr,vol,sep)
yctr= 0;
% vol = 4/3*pi*r^3
r1 = (3/4*vol/pi)^(1/3);
r2 = (3/4*vol/2/pi)^(1/3);
balltxt=   'solid ball%d   = sphere(%f,%f,%f;%f);\n';
ball1 = sprintf(balltxt,1,xctr    ,yctr,zctr,r1); 
ball2 = sprintf(balltxt,2,xctr+sep,yctr,zctr,r2); 
ball3 = sprintf(balltxt,3,xctr-sep,yctr,zctr,r2); 

shape_str = [ ...
   'solid top    = plane(0,0,0;0,0,1);\n' ...
   ball1, ball2, ball3, ...
   'solid ball21 = ball2 and      ball1 ; tlo ball21;\n', ...
   'solid ball2_ = ball2 and (not ball1); tlo ball2_;\n', ...
   'solid ball31 = ball3 and      ball1 ; tlo ball31;\n', ...
   'solid ball3_ = ball3 and (not ball1); tlo ball3_;\n', ...
   'solid ball1_ = ball1 and (not ball2) and (not ball3); tlo ball1_;\n', ...
   'solid balls = ball1 or ball2 or ball3;\n', ...
   'solid mainobj= top and orthobrick(-3,-3,-4;3,3,0) and not balls -maxh=1.0;\n'];
[elec_pos_x,elec_pos_y] = meshgrid(linspace( -1.5,1.5,5),linspace(-2,2,7));
elec_pos = [  elec_pos_x(:), elec_pos_y(:), ones(size(elec_pos_x(:)))*[0,0,0,1] ];
elec_shape=[0.1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

mi1 = [fmdl.mat_idx{1}; fmdl.mat_idx{3}; fmdl.mat_idx{5}];
mi2 = [fmdl.mat_idx{1}; fmdl.mat_idx{2};  ...
       fmdl.mat_idx{3}; fmdl.mat_idx{4}]; 
return

img.elem_data(fmdl.mat_idx{1}) = 1.1;
img.elem_data(fmdl.mat_idx{2}) = 0.9;
img.elem_data(fmdl.mat_idx{3}) = 1.1;
img.elem_data(fmdl.mat_idx{4}) = 0.9;
show_fem(img);
