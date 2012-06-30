% 3D Simulatulation example of system evaluation
% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3
function [DET] = sim_Det_3D

in_plane_Def=0; % in plane or off-plane
Nelec= 16; Dims = 3; nrings = 1;
opt = {'no_meas_current'};
sp = 1; mp= 1;  % stimulation, measurement

[img, vol, boxidx] = calc_fwd_cube( Dims, Nelec,in_plane_Def);
stim = mk_stim_patterns(Nelec, nrings, [0,sp], [0,mp], opt, 1);
DET = distinguishability_calc( img, stim, boxidx, vol);
plott(DET);
end

function SNR = distinguishability_calc( img, stim, boxidx, vol)
   img.fwd_model.stimulation = stim;
   J  = calc_jacobian( img );
   for ri = 1:length(boxidx);
      idx = boxidx{ri};
      Jr = J*idx;
      SNR(ri) = 1e6*sqrt( Jr'*Jr );
   end
 end

function [img, vol, boxidx] = calc_fwd_cube(Dims, Nelec,in_plane_Def)
% cube moved to 7 position
if in_plane_Def  %% in plane

extra={'boxes', ...
    ['solid box1= orthobrick(-0.1,-0.85,0.9; 0.1,-0.65,1.1) -maxh=0.05;' ...
    'solid box2= orthobrick(-0.1, -0.6,0.9; 0.1, -0.4,1.1) -maxh=0.05;' ...
    'solid box3= orthobrick(-0.1, -0.35,0.9; 0.1,-0.15,1.1) -maxh=0.05;' ...
    'solid box4= orthobrick(-0.1,-0.1,0.9; 0.1,0.1,1.1) -maxh=0.05;' ...   
    'solid box5= orthobrick(-0.1, 0.15,0.9; 0.1,0.35,1.1) -maxh=0.05;' ...
    'solid box6= orthobrick(-0.1, 0.4,0.9; 0.1,0.6,1.1) -maxh=0.05;' ...
    'solid box7= orthobrick(-0.1, 0.65,0.9; 0.1,0.85,1.1) -maxh=0.05;' ...
 'solid boxes = box1 or box2 or box3 or box4 or box5 or box6 or box7;']};
else  %% off plane
extra={'boxes', ...
    ['solid box1= orthobrick(-0.1,-0.85,1.4; 0.1,-0.65,1.6) -maxh=0.05;' ...
    'solid box2= orthobrick(-0.1, -0.6,1.4; 0.1, -0.4,1.6) -maxh=0.05;' ...
    'solid box3= orthobrick(-0.1, -0.35,1.4; 0.1,-0.15,1.6) -maxh=0.05;' ...
    'solid box4= orthobrick(-0.1,-0.1,1.4; 0.1,0.1,1.6) -maxh=0.05;' ...   
    'solid box5= orthobrick(-0.1, 0.15,1.4; 0.1,0.35,1.6) -maxh=0.05;' ...
    'solid box6= orthobrick(-0.1, 0.4,1.4; 0.1,0.6,1.6) -maxh=0.05;' ...
    'solid box7= orthobrick(-0.1, 0.65,1.4; 0.1,0.85,1.6) -maxh=0.05;' ...
 'solid boxes = box1 or box2 or box3 or box4 or box5 or box6 or box7;']};
end

switch Dims;  % dims = 2 (2D) or 3 (3d)
    case 3;             
        fmdl= ng_mk_cyl_models(2,[Nelec,1],[0.1,0,0.03],extra);
 %show_fem(fmdl, [1 1 0])
    case 2;
        fmdl= ng_mk_cyl_models(0,Nelec,[0.1,0,0.01],extra);
    otherwise; error('huh?');
end

img= mk_image(fmdl,1);
idx = fmdl.mat_idx{2};
vol= get_elem_volume(fmdl);
ctrs = interp_mesh(fmdl); ctrs = ctrs( idx,: );
rctrs = sqrt(ctrs(:,1).^2 + ctrs(:,2).^2) ;
x= ctrs(:,1); y= ctrs(:,2); z= ctrs(:,3);
contr = 5*1e-3;
area0 = sum(vol(idx));

indx = zeros(size(vol)); list = idx(y < -0.625);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{1} = indx.*vol;

% img.elem_data = boxidx{1};
% show_slices(img,1)

indx = zeros(size(vol)); list = idx( y < -0.375 & y > -0.625);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{2} = indx.*vol;

indx = zeros(size(vol)); list = idx( y < -0.125 & y > -0.375);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{3} = indx.*vol;

indx = zeros(size(vol)); list = idx(y > -0.125 & y < 0.125);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{4} = indx.*vol;

indx = zeros(size(vol)); list = idx( y > 0.125 & y < 0.375);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{5} = indx.*vol;

indx = zeros(size(vol)); list = idx( y > 0.375 & y < 0.625);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{6} = indx.*vol;

indx = zeros(size(vol)); list = idx( y > 0.625);
area1 = sum(vol(list));
indx(list)=+ contr*area0/area1;
boxidx{7} = indx.*vol;
end

function plott(DET)
ax = [-0.75, -0.5, -0.2, 0, 0.2, 0.5, 0.75];
plot(ax, DET/max(DET),'linewidth',4)
ylabel('Normalized Detectability {\it z}');
xlabel('Radius');
end
