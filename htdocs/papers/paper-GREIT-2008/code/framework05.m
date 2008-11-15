% Test RM $Id$

% EXCLUDE MEASURES AT ELECTRODES
[x,y]= meshgrid(1:16,1:16);
idx= abs(x-y)>1 & abs(x-y)<15;

% LOAD SOME TEST DATA
load iirc_data_2006
v_reference= - real(v_reference(idx,:));
v_rotate   = - real(v_rotate(idx,:));

load RM_framework_example RM map;
for k=1:4;
   dv = v_rotate(:,k*2) - v_reference;
   img = reshape( RM*dv, 32,32); % reconstruction
   img(~map) = NaN;              % background
   imgs(:,:,k) = img;
end

imagesc(reshape(imgs,32,4*32))
axis off; axis equal; colormap(gray(256))
print -dpng -r100 framework05a.png
