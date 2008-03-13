% Lung images
% $Id: tutorial410c.m,v 1.2 2008-03-13 22:08:30 aadler Exp $

raster_img= calc_slices(img);
raster_img(isnan(raster_img))=0;
% define roi as whole image
s_ri = size(raster_img);
roi = ones(s_ri(1:2));

for i= 1:s_ri(3)
  sig(i)= sum(sum(raster_img(:,:,i) .* roi));
end

subplot(221)
plot( ((1:s_ri(3))-1)*5, sig/sig(1))
xlabel('minutes after drink')
ylabel('normalized conductivity');

print -r125 -dpng tutorial410c.png;
