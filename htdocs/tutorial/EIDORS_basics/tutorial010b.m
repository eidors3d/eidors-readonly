% Simulate EIT data
% $Id: tutorial010b.m,v 1.1 2007-06-15 18:17:51 aadler Exp $

sim_img= eidors_obj('image', 'stimulation image');
sim_img.fwd_model= imdl_3d.fwd_model;

% set homogeneous conductivity and simulate
sim_img.elem_data= ones( size(sim_img.fwd_model.elems,1) ,1);
homg_data=fwd_solve( sim_img );

% set inhomogeneous conductivity and simulate
sim_img.elem_data([390,391,393,396,402,478,479,480,484,486, ...
                   664,665,666,667,668,670,671,672,676,677, ...
                   678,755,760,761])= 1.15;
sim_img.elem_data([318,319,321,324,330,439,440,441,445,447, ...
                   592,593,594,595,596,598,599,600,604,605, ...
                   606,716,721,722])= 0.8;
inh_data=fwd_solve( sim_img );

clf;subplot(211);

xax= 1:length(homg_data.meas);
hh= plotyy(xax,[homg_data.meas, inh_data.meas], ...
           xax, homg_data.meas- inh_data.meas );

set(hh,'Xlim',[1,max(xax)]);
print -r75 -dpng tutorial010b.png
