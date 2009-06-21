function igt = img2igt(img)
%IMG2IGT returns an IGT-compatible image matrix from any EIDORS
% IMG struct. 
%
% IGT = IMG2IGT(IMG) returns a vector NFrames-by-912.
%
% WARNING: When the mesh stored in fwd_model is not a resterised 32-by-32
% matrix, this conversion results in a loss of quality.
% 
% See also EIDORS_SAVEIMG.

% (C) 2009 by Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$ 

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

% create IGT vector
igt=zeros(size(img.elem_data,2),912);
for i = 1:32
    for j = 1:32
        if igt_m(i,j)>0 
            igt(:,igt_m(i,j)) = greit_m(i,j,:);
        end
    end
end
