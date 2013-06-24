
SimModel=3;
LungFlag=0; % if lung areas should be simulated
switch SimModel
    case 1
        %***(1) Cylindrical,
        % extra={'ball','solid ball = sphere(0,0.5,1;0.2) or sphere(0,-0.5,1;0.2);'};
        % extra={'ball','solid ball = cylinder(0.2,-0.2,0;0.2,-0.2,1;0.2) and orthobrick(-1,-1,0;1,1,1);'};
        if LungFlag==1
            extra={'lung','solid lung = ellipsoid(0.4,0,1; 0.25,0,0; 0,0.6,0; 0,0,0.6) or ellipsoid(-0.4,0,1; 0.25,0,0; 0,0.6,0; 0,0,0.6);'}; % x,y,z; l/2,0,0; 0,h/2,0; 0,0,w/2
            fmdl= ng_mk_cyl_models([2,1,.15],[16,1.0],[0.1],extra);
        else
%             fmdl= ng_mk_cyl_models([0.3,1,.05],[16,0.15],[0.1]);  % small Z value better figures of merit
                fmdl= ng_mk_cyl_models([2,1,.15],[16,1.0],[0.1]);
        end
        show_fem(fmdl);
    case 2
        %***(2) elliptical,
        % Put two ellipsoid into the elliptical cylinder
        extra={'lung','solid lung = ellipsoid(0.4,0,1; 0.25,0,0; 0,0.45,0; 0,0,0.6) or ellipsoid(-0.4,0,1; 0.25,0,0; 0,0.45,0; 0,0,0.6);'}; % x,y,z; l/2,0,0; 0,h/2,0; 0,0,w/2
        fmdl= ng_mk_ellip_models([2,1,0.7],[16,1],[0.1],extra);
        show_fem(fmdl);
    case 3
        %***(3) thoracic (pig/human) forward models
        % get contours
        thorax = shape_library('get','adult_male','boundary');
        rlung  = shape_library('get','adult_male','right_lung');
        llung  = shape_library('get','adult_male','left_lung');
        % one could also run:
        % shape_library('get','adult_male');
        % to get all the info at once in a struct
        if LungFlag==1
            shape = { 2,                      % height
                {thorax, rlung, llung}, % contours
                [4,50],                 % perform smoothing with 50 points
                0.15};                  % large maxh (coarse mesh)
        else
            shape = { 2,                      % height
                {thorax}, % contours
                [4,50],                 % perform smoothing with 50 points
                0.15};                  % large maxh (coarse mesh)
        end
        elec_pos = [ 16,                  % number of elecs per plane
            1,                   % equidistant spacing
            1]';               % a single z-plane
        
        elec_shape = [0.1,               % radius
            0,                  % circular electrode
            0.05 ]';             % maxh (electrode refinement)
        
        fmdl = ng_mk_extruded_model(shape, elec_pos, elec_shape);
        % this similar model is also available as:
        % fmdl = mk_library_model('adult_male_16el_lungs');
        
        
%         img=mk_image(fmdl,1);
        % img.elem_data(fmdl.mat_idx{2})= 0.3; % rlung
        % img.elem_data(fmdl.mat_idx{3})= 0.3; % llung
        
%         clf; show_fem(img); view(0,70);


    case 4
%%(4) maybe add the water tank phantom

end
%%
stim = mk_stim_patterns(16,1,[0,1],[0,1],{}); 
fmdl.stimulation = stim;
img = mk_image(fmdl, 1); 
vh = fwd_solve(img);
if LungFlag==1
    if SimModel==3
        img.elem_data(fmdl.mat_idx{2})= 0.3; % rlung
        img.elem_data(fmdl.mat_idx{3})= 0.3; % llung
    else
        img.elem_data(fmdl.mat_idx{2})= 2.0;
    end
end
vi = fwd_solve(img);
show_fem(img);

%%
fmdl.normalize_measurements = 1;

opt.imgsz = [32 32];
opt.square_pixels = 0;
opt.Nsim = 1000;
opt.distr = 3;
opt.target_size=0.05;
% opt.target_plane
% opt.target_offset
opt.noise_figure = 0.5;

i_gr = mk_GREIT_model(img, 0.25, [], opt);

Test = 2; % different simulation patern
if Test==1
    %% Test 1
    % [r, params] =  test_performance( i_gr );
    Nsim = 100;
%     l =  linspace(0,1,Nsim);
    l =  linspace(0,0.9,Nsim);
    ctr = [0,0, mean(i_gr.fwd_model.nodes(:,3))];  % Assume x,y centre is zero
    maxx = max(abs(i_gr.fwd_model.nodes(:,1) - ctr(1)));
    maxy = max(abs(i_gr.fwd_model.nodes(:,2) - ctr(2)));
    r = 0.05; % target radius
    th=  l*4321; % want object to jump around in radius
    Zpos=ctr(3)*ones(1,Nsim);
%     xyzr = [maxx*l.*cos(th); maxy*l.*sin(th); ones(1, Nsim); 0.05*ones(1, Nsim)];
    xyzr = [maxx*l.*cos(th); maxy*l.*sin(th); Zpos; r*mean([maxx,maxy])*ones(1, Nsim)];
    
elseif Test==2
    %% Test 2 for thorax SimModel=3;
    r = 0.05; % target radius
    Npos = 50; % number of positions
    Xpos =  linspace(-0.9,0.9,Npos); % positions to simulated along x-axis
    Xpos2 = -0.2*ones(1,Npos); % along y-axis
    Xpos=[Xpos Xpos Xpos2]; % 1 along x; 2 along x, off plane; 3 along y
    Ypos = zeros(1,Npos);
    Ypos2 = linspace(-0.55,0.6,Npos);
    Ypos=[Ypos Ypos Ypos2];
    Zpos = ones(1,Npos);  %% for off-plane, adjust the level (*1.5)
    Zpos = [Zpos Zpos*0.7 Zpos];
%     xyzr = [Xpos; Ypos; 0.15*Zpos; r*ones(1,Npos*3)];  % corelate to small Z in fmdl
    xyzr = [Xpos; Ypos; Zpos; r*ones(1,Npos*3)];
    [xyzr]=del_out_map(img,64,xyzr);
else
    %% Test 3 old Test_performance from Bartek
    N = 16;
    bnd_nodes = unique(i_gr.fwd_model.boundary);
    min_bb = min(i_gr.fwd_model.nodes(bnd_nodes,:));
    max_bb = max(i_gr.fwd_model.nodes(bnd_nodes,:));
    h = mean([min_bb(3),max_bb(3)]);            % better way to calculate Zpos
    r = 0.025 * max(max_bb(1:2) - min_bb(1:2));
    xspace = linspace(min_bb(1),max_bb(1),N);
    yspace = linspace(min_bb(2),max_bb(2),N);
    [X Y] = meshgrid(xspace,yspace);
    img.calc_colours.npoints = N;
    M = calc_slices(img,1);
    IN = M==1;
    n_px = nnz(IN);
    OUT = isnan(M);
    Zpos=h*ones(1,n_px);
%     xyzr = [X(IN)'; Y(IN)'; Zpos; r*ones(1,n_px); ];
    xyzr = [0.8*X(IN)'; 0.8*Y(IN)'; Zpos; r*ones(1,n_px); ];
end

[vh,vi] = simulate_movement(img, xyzr);
imgr= inv_solve(i_gr,vh,vi);

% calculate the GREIT parameters
figure
levels =[inf,inf,Zpos];
imgr.show_slices.img_cols =10;
show_slices(imgr, levels);
imgr.calc_colours.npoints = 64;
imgr.calc_slices.levels=levels;
params = eval_GREIT_fig_merit(imgr, xyzr);
figure
p_names = {'AR','PE','RES','SD','RNG'};
for i=1:5; subplot(5,1,i);
    plot(params(i,:)); ylabel(p_names{i});
    mean_save(i)=mean(abs(params(i,:)));
    std_save(i)=std(abs(params(i,:)));
    
end

