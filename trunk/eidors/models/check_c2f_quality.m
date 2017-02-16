function [ f_frac, c_frac ] = check_c2f_quality( f_mdl, c_mdl, c2f )
%CHECK_C2F_QUALITY Check the coarse2fine mapping between two models
%   CHECK_C2F_QUALITY(F_MDL, C_MDL, C2F) creates a plot for assesing the
%   quality of a coars2fine mapping.
%   A coarse2fine mapping matrix C2F between a fine model F_MDL and a 
%   coarse model C_MDL contains in each element C2F(i,j) the fraction of 
%   F_MDL element i contained in C_MDL element j. CHECK_C2F_QUALITY 
%   calculates the fraction of the volume of each element of both models
%   that is covered by the C2F mapping. Ideally, the result is always 1,
%   other than in places where models don't overlap.
%
%   [ f_frac, c_frac ] = CHECK_C2F_QUALITY(...) returns the calculated
%   volumes for further analysis.

% (C) 2017 Bartlomiej Grychtol 
% License: GPL version 2 or 3
% $Id:$

c_vol = get_elem_volume(c_mdl);
f_vol = get_elem_volume(f_mdl);

% C2F(i,j) is the fraction of f_mdl element i contained in c_mdl element j

% Fraction of f_mdl covered by c2f:
f_frac = full(sum(c2f,2));
f_img = mk_image(f_mdl, f_frac);

% Fraction of c_mdl elements covered by c2f
c_frac = (c2f'*f_vol) ./ c_vol;
c_img = mk_image(c_mdl, c_frac);

subplot(221)
show_fem(f_img)
title('Fine model coverage')
% hold on
% h = show_fem(c_mdl);
% set(h,'EdgeColor','b')
eidors_colourbar(f_img);
hold off

subplot(222)
hist(f_frac,1000);
title('Fine model coverage')

subplot(223)
show_fem(c_img)
% hold on
% h = show_fem(f_mdl);
% set(h,'EdgeColor','b')
eidors_colourbar(c_img);
hold off
title('Coarse model coverage')


subplot(224)
hist(c_frac,1000);
title('Coarse model coverage');
end

