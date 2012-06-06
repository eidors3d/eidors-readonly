function hm_img = calc_hm_set(img,frac)
% hm_img= CALC_HA_SET(img)
% hm_img is szxszxNimg. It is 1 inside the Half Ampl Set
% frac is the fraction of maximum (0.5 or 0.25)
% hm_img expects conductive changes. Use calc_hm_set(-img,frac) for non-c

% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

sz = size(img,1);
[x,y]=meshgrid(linspace(-1,1,sz),linspace(-1,1,sz)); map = x.^2+y.^2<1.1;

hm_img = logical(zeros(size(img)));
for i=1:size(img,3);
   imi = img(:,:,i); imi= imi(map);

   hmi= logical(zeros(sz));
   hmi(map) = imi >= (max(imi) * frac);
   hm_img(:,:,i) = hmi;
end
