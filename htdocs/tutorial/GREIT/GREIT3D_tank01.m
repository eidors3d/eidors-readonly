%STEP 1: Simple model with vertical targets
posns= linspace(0.5,3.5,13);
str=''; for i=1:length(posns);
   extra{i} = sprintf('ball%03d',round(posns(i)*100));
   str = [str,sprintf('solid %s=sphere(0.5,0,%f;0.1); ', extra{i}, posns(i))];
end
extra{i+1} = str;
fmdl= ng_mk_cyl_models([4,1,.2],[16,1.5,2.5],[0.05],extra); 
fmdl = mdl_normalize(fmdl, 0);
[~,fmdl] = elec_rearrange([16,2],'square', fmdl);

img= mk_image(fmdl,1);
img.elem_data(vertcat(fmdl.mat_idx{2:end}))= 2;
opt.viewpoint = struct('az',0,'el',10);
clf;show_fem_enhanced(img,opt);

print_convert GREIT3D_tank01a.jpg
