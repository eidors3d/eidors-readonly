function ha_img = calc_ha_set(img)
% ha_img= CALC_HA_SET(img)
% ha_img is 32x32xNimg. It is 1 inside the Half Ampl Set

% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

[x,y]=meshgrid(linspace(-1,1,32),linspace(-1,1,32));
map = x.^2+y.^2<1.1;


ha_img = logical(zeros(size(img)));
for i=1:size(img,3);
   imi = img(:,:,i); imi= imi(map);
   ss_imi= sum(imi);
   [s_imi, idx] = sort(imi*sign(ss_imi));
   [jnk, ff] = min( abs( cumsum(s_imi) - abs(ss_imi)/4 ) );

   hai = logical(zeros(sum(map(:)),1));
   hai(idx(ff(1)+1:end))= 1;
   hai2= logical(zeros(32));
   hai2(map) = hai;
   ha_img(:,:,i) = hai2;
end
