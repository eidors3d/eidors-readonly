function igt = img2igt(img)
% IMG2IGT - returns an IGT-compatible measurement matrix from any EIDORS
% IMG struct.

% (C) 2009 by Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$ 

% loc.name = 'Dummy vector 1:856';
% loc.elem_data = 1:856;
% loc.fwd_model = img.fwd_model; 
% loc.type = 'image';

img.calc_colours.npoints=32;
greit_m = calc_slices(img);
greit_m(isnan(greit_m)) = 0;

ind=[8:25,39:58,70:91,101:124,132:157,163:190,194:223,225:800, ...
        802:831,835:862,868:893,901:924,934:955,967:986,1000:1017];
igt_m = zeros(1,1024);
igt_m(ind) = 1:912;
igt_m = reshape(igt_m,32,32)';

igt=zeros(size(img.elem_data,2),912);

for i = 1:32
    for j = 1:32
        if igt_m(i,j)>0 %&& greit_m(i,j)>0
            %igt(:,igt_m(i,j)) = img.elem_data(greit_m(i,j),:);
            igt(:,igt_m(i,j)) = greit_m(i,j,:);
        end
    end
end
