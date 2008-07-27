function hm_img = calc_hm_set(img,frac)
% hm_img= CALC_HA_SET(img)
% hm_img is 32x32xNimg. It is 1 inside the Half Ampl Set
% frac is the fraction of maximum (0.5 or 0.25)

% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

[x,y]=meshgrid(linspace(-1,1,32),linspace(-1,1,32)); map = x.^2+y.^2<1.1;

hm_img = logical(zeros(size(img)));
for i=1:size(img,3);
   imi = img(:,:,i); imi= imi(map);
   if abs(max(imi)) < abs(min(imi)); imi=-imi; end

   hmi= logical(zeros(32));
   hmi(map) = imi >= max(imi) * frac;
   hm_img(:,:,i) = hmi;
end
