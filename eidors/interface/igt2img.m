function img = igt2img(igt)
%IGT2IMG constructs an EIDORS IMG struct from an IGT frames-by-912 matrix. 
% IMG = IGT2IMG(IGT)
% 
% See also IMG2IGT, EIDORS_READIMG, EIDORS_SAVEIMG.

% (C) 2009 by Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$


img.name = 'Read in from an igt file';
img.type = 'image';
img.elem_data = 1:912;
tempmdl = mk_common_gridmdl('backproj');
img.fwd_model = tempmdl.fwd_model;

% render the image on a 32-by-32 grid
img.calc_colours.npoints=32;
greit_m = calc_slices(img);
greit_m(isnan(greit_m)) = 0;

% create an IGT mask
ind=[8:25,39:58,70:91,101:124,132:157,163:190,194:223,225:800, ...
        802:831,835:862,868:893,901:924,934:955,967:986,1000:1017];
igt_m = zeros(1,1024);
igt_m(ind) = 1:912;
igt_m = reshape(igt_m,32,32)';

% create data vector
data=zeros(size(igt,1),912);
for i = 1:32
    for j = 1:32
        if igt_m(i,j)>0 
            data(:,greit_m(i,j)) = igt(:,igt_m(i,j));
        end
    end
end

img.elem_data = data;
