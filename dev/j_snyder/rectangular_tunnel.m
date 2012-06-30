x_distance = 2.5;
inner_rad = 1;
outer_rad = 5;
y_excursion = 2.5;

targposes = [];
estposes = [];

N_elec = 16;

shape_strf = ['solid incyl  = cylinder (0,0,-1; 0,0,0; ' num2str(inner_rad) ') -maxh=0.2; \n' ...
             'solid farblock = orthobrick (-' num2str(outer_rad) ',-' num2str(outer_rad) ',-1.25; ' num2str(outer_rad) ',' num2str(outer_rad) ',1.25 ) -maxh=0.5; \n' ...
             'solid pl1    =  plane(0,0,-1.25;0,0,-1);\n' ...
             'solid pl2    =  plane(0,0,1.25; 0,0,1);\n' ...
             'solid mainobj= pl1 and pl2 and farblock and not incyl; \n'];
fmdl = eidors_obj('get-cache', shape_strf, 'forward_model',N_elec);
if isempty(fmdl)
    th= linspace(0,2*pi,N_elec+1)'; th(end)=[];
    cth= cos(th); sth=sin(th); zth= zeros(size(th));
    elec_pos = [cth, sth, zth, cth, sth, zth];
    elec_shape=[0.01];
    elec_obj = 'incyl';
    [fmdl, mat_idx] = ng_mk_gen_models(shape_strf, elec_pos, elec_shape, elec_obj);

    eidors_obj('set-cache', shape_strf, 'forward_model', fmdl, N_elec);
end


% Create a 2D Slice by making a flattened version of the 3D model
shape_strc = ['solid incyl  = cylinder (0,0,-1; 0,0,0; ' num2str(inner_rad) ') -maxh=0.1; \n' ...
             'solid farblock = orthobrick (-' num2str(outer_rad) ',-' num2str(outer_rad) ',-0.1; ' num2str(outer_rad) ',' num2str(outer_rad) ',0 ) -maxh=0.25; \n' ...
             'solid pl1    =  plane(0,0,-0.1;0,0,-1);\n' ...
             'solid pl2    =  plane(0,0,0; 0,0,1);\n' ...
             'solid mainobj= pl1 and pl2 and farblock and not incyl; \n'];
  
cmdl = eidors_obj('get-cache', shape_strc, 'reconst_model',N_elec);
if isempty(cmdl)
    th= linspace(0,2*pi,N_elec+1)'; th(end)=[];
    cth= cos(th); sth=sin(th); zth= zeros(size(th))-0.05;
    elec_pos = [cth, sth, zth, cth, sth, zth];
    elec_shape=[0.03 0.2 0.1];
    elec_obj = 'incyl';
    [f2mdl, mat_idx] = ng_mk_gen_models(shape_strc, elec_pos, elec_shape, elec_obj);
    cmdl = mdl2d_from3d(f2mdl,mat_idx);
    
    eidors_obj('set-cache', shape_strc, 'reconst_model', fmdl, N_elec);
end

% Simulation protocol. TODO: we need a geophysical stim protocol
stim = mk_stim_patterns(N_elec, 1, [0,5], [0,5], {'no_meas_current'},1);
fmdl.stimulation = stim;
cond_mdl = .02; % in S/m units
img = mk_image( fmdl, cond_mdl); 
vs_h = fwd_solve( img);

i = 1;
for ypos = -y_excursion:0.2:y_excursion
img.elem_data = cond_mdl*(1 + 1000*mk_c2f_circ_mapping(fmdl, [x_distance;ypos;0;0.2]) );
vs_i = fwd_solve( img);

%vs_i= eidors_obj('data', 'noisy data', ...
%                           'meas', vs_i.meas + 1e-3*randn(size(vs_i.meas,1),1));
figure(1); clf;
show_fem(img);
view(0,90);
 



imdl = select_imdl( fmdl, {'Basic GN dif'});
imdl.rec_model = cmdl;
imdl.jacobian_bkgnd.value = cond_mdl;

% Do coarse2fine mapping. Rotate mdl to z dirn
f1mdl = fmdl;
f1mdl.mk_coarse_fine_mapping.z_depth = 1;
c2f= mk_coarse_fine_mapping( f1mdl, cmdl);
imdl.fwd_model.coarse2fine = c2f;


imdl.solve=    @aa_inv_solve;
imdl.hyperparameter.value = 0.5e-2;

imgr = inv_solve( imdl, vs_h, vs_i );

imgr.calc_colours.npoints= 512;
[val idx] = max(imgr.elem_data);
xyest = mean(imgr.fwd_model.nodes(imgr.fwd_model.elems(idx,:),:));
estposes = [estposes; xyest];
figure(2); clf;
show_fem(imgr);
hold on;
targposes = [targposes;x_distance ypos];

plot(targposes(end,1),targposes(end,2),'kx','MarkerSize',20,'LineWidth',2)
plot(estposes(end,1),estposes(end,2),'ko','MarkerSize',20,'LineWidth',2)
print('-dpng',sprintf('tunnelsim_frame_%08d',i))

i = i + 1;
end
