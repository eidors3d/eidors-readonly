
SimModel=3;
LungFlag=0; % if lung areas should be simulated
switch SimModel
    case 1
        %***(1) Cylindrical,
        % extra={'ball','solid ball = sphere(0,0.5,1;0.2) or sphere(0,-0.5,1;0.2);'};
%         extra={'ball','solid ball = cylinder(0.2,-0.2,0;0.2,-0.2,1;0.2) and orthobrick(-1,-1,0;1,1,1);'};
        if LungFlag==1
%             extra={'lung','solid lung = ellipsoid(0.4,0,1; 0.25,0,0;
%             0,0.6,0; 0,0,0.6) or ellipsoid(-0.4,0,1; 0.25,0,0; 0,0.6,0; 0,0,0.6);'}; % x,y,z; l/2,0,0; 0,h/2,0; 0,0,w/2
            
% ******   thorax +0.25
            extra={'ball1','ball2',['solid ball1 = ellipsoid(0.47,0,1.5; 0.33,0,0; 0,0.6,0; 0,0,0.6);' 'solid ball2 = ellipsoid(-0.47,0,1.5; 0.4,0,0; 0,0.6,0; 0,0,0.6);']}; % x,y,z; l/2,0,0; 0,h/2,0; 0,0,w/2
            fmdl= ng_mk_cyl_models([3,1.25,.15],[16,1.5],[0.1],extra);
            
        else
%             fmdl= ng_mk_cyl_models([0.3,1,.05],[16,0.15],[0.1]);  % small Z value better figures of merit
                fmdl= ng_mk_cyl_models([3,1,.15],[16,1.5],[0.1]);
        end
        show_fem(fmdl);
    case 2
        %***(2) elliptical,
        % Put two ellipsoid into the elliptical cylinder
        %symetric
        extra={'lung1','solid lung1 = ellipsoid(0.4,0,1; 0.25,0,0; 0,0.45,0; 0,0,0.6) or ellipsoid(-0.4,0,1; 0.25,0,0; 0,0.45,0; 0,0,0.6);'}; % x,y,z; l/2,0,0; 0,h/2,0; 0,0,w/2  symetric
        % asymetric
%         extra={'ball1','ball2',['solid ob = orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;' 'solid ball1 = ellipsoid(0.4,0,1; 0.25,0,0; 0,0.45,0; 0,0,0.6) and ob;' 'solid ball2 = ellipsoid(-0.4,0,1; 0.2,0,0; 0,0.45,0; 0,0,0.65) and ob;']}; % x,y,z; l/2,0,0; 0,h/2,0; 0,0,w/2

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
            shape = { 3,                      % height
                {thorax, llung, rlung}, % contours
                [4,50],                 % perform smoothing with 50 points
                0.15};                  % large maxh (coarse mesh)
        else
            shape = { 3,                      % height
                {thorax}, % contours
                [4,50],                 % perform smoothing with 50 points
                0.15};                  % large maxh (coarse mesh)
        end
        elec_pos = [ 16,                  % number of elecs per plane
            1,                   % equidistant spacing
            1.5]';               % a single z-plane
        
        elec_shape = [0.1,               % radius
            0,                  % circular electrode
            0.05 ]';             % maxh (electrode refinement)
        
        fmdl = ng_mk_extruded_model(shape, elec_pos, elec_shape);
       
        show_fem(fmdl);

    case 4
%%(4) maybe add the water tank phantom

end
%%
stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current','rotate_meas'}); 
% stim = mk_stim_patterns(16,1,[0,1],[0,1],{}); 
fmdl.stimulation = stim;
fmdl.normalize_measurements = 1;
img_h = mk_image(fmdl, 1); 
img=img_h;
vh = fwd_solve(img);
if LungFlag==1
    if SimModel==3 || SimModel==1
        img.elem_data(fmdl.mat_idx{2})= 0.6; % llung
        img.elem_data(fmdl.mat_idx{3})= 0.6; % rlung
    else
        img.elem_data(fmdl.mat_idx{2})= 1.4;
    end
end
vi = fwd_solve(img);
show_fem(img);

%%
opt.imgsz = [32 32];
opt.square_pixels = 0;
opt.Nsim = 2000;
opt.distr = 6;
opt.target_size=0.05;
% opt.target_plane
% opt.target_offset
opt.noise_figure = 0.5;

i_gr = mk_GREIT_model(fmdl, 0.25, [], opt);

%%
Test = 6; % different simulation patern
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
    Zpos = ones(1,Npos);  %% for off-plane, adjust the level (*0.7)
    Zpos = [Zpos Zpos*0.7 Zpos];
%     xyzr = [Xpos; Ypos; 0.15*Zpos; r*ones(1,Npos*3)];  % corelate to small Z in fmdl
    xyzr = [Xpos; Ypos; Zpos; r*ones(1,Npos*3)];
    [xyzr]=del_out_map(img,64,xyzr);
elseif Test==3
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
%     [xyzr]=del_out_map(img,64,xyzr);
    xyzr = [[0;0.9*X(IN)]'; [0;0.9*Y(IN)]'; [h,Zpos]; r*ones(1,n_px+1); ]; % add zero point for reference
elseif Test == 4
    %% Test 4 old Test_performance from Bartek + off plane objects
    N = 32;
    bnd_nodes = unique(i_gr.fwd_model.boundary);
    min_bb = min(i_gr.fwd_model.nodes(bnd_nodes,:));
    max_bb = max(i_gr.fwd_model.nodes(bnd_nodes,:));
    h = mean([min_bb(3),max_bb(3)]);            % better way to calculate Zpos
    r = 0.025 * max(max_bb(1:2) - min_bb(1:2));
    xspace = linspace(min_bb(1),max_bb(1),N);
    yspace = linspace(min_bb(2),max_bb(2),N);
    [X Y] = meshgrid(xspace,yspace);
    img1 = mk_image(i_gr.fwd_model,1);
    img1.calc_colours.npoints = N;
    M = calc_slices(img1,1);
    IN = M==1;
    n_px = nnz(IN);
    OUT = isnan(M);
%     Zpos=h*ones(1,n_px);
    %     xyzr = [X(IN)'; Y(IN)'; Zpos; r*ones(1,n_px); ];
    Z_layers=1;
    temp_X=0.9*X(IN);
    temp_Y=0.9*Y(IN);

    X_all=repmat(temp_X,Z_layers,1);
    Y_all=repmat(temp_Y,Z_layers,1);
    if Z_layers~=1
        hspace=linspace(h, max_bb(3),Z_layers);
        Z_all=hspace'*ones(1,n_px);
        Z_all=reshape(Z_all',size(Z_all,1)*size(Z_all,2),1);
    else
        Z_all=h*ones(1,n_px)';
    end
    %     xyzr = [X_all'; Y_all'; Z_all'; r*ones(1,length(Z_all)); ];
    xyzr = [[0;X_all]'; [0;Y_all]'; [h;Z_all]'; r*ones(1,length(Z_all)+1); ]; %add a central point. assume x=0,y=0;
    Zpos=Z_all';

elseif Test == 5
    N = 32;
     % this bit is copied from mdl_slice_mapper to make sure
     xmin = min(i_gr.fwd_model.nodes(:,1));    xmax = max(i_gr.fwd_model.nodes(:,1));
     xmean= mean([xmin,xmax]); xrange= xmax-xmin;

     ymin = min(i_gr.fwd_model.nodes(:,2));    ymax = max(i_gr.fwd_model.nodes(:,2));
     ymean= mean([ymin,ymax]); yrange= ymax-ymin;

     range= max([xrange, yrange]);
     xspace = linspace( xmean - range*0.5, xmean + range*0.5, N );
     yspace = linspace( ymean + range*0.5, ymean - range*0.5, N );
     % end of copied block
    bnd_nodes = unique(i_gr.fwd_model.boundary);
    min_bb = min(i_gr.fwd_model.nodes(bnd_nodes,:));
    max_bb = max(i_gr.fwd_model.nodes(bnd_nodes,:));
    h = mean([min_bb(3),max_bb(3)]);
    r = 0.025 * range;
%     xspace = linspace(min_bb(1),max_bb(1),N);
%     yspace = linspace(max_bb(2),min_bb(2),N);
    [X Y] = meshgrid(xspace,yspace);
    img1 = mk_image(i_gr.fwd_model,1);
    img1.calc_colours.npoints = N;
    M = calc_slices(img1,1);
    IN = M==1;
    n_px = nnz(IN);
    OUT = isnan(M);
    Zpos = h*ones(1,n_px);
    xyzr = [0.9*X(IN)'; 0.9*Y(IN)'; Zpos;  r*ones(1,n_px); ];
    
elseif Test == 6
    %% Test 6 objects first from inside then outside + off plane objects
    N = 32;
    bnd_nodes = unique(i_gr.fwd_model.boundary);
    min_bb = min(i_gr.fwd_model.nodes(bnd_nodes,:));
    max_bb = max(i_gr.fwd_model.nodes(bnd_nodes,:));
    h = mean([min_bb(3),max_bb(3)]);            % better way to calculate Zpos
    r = 0.025 * max(max_bb(1:2) - min_bb(1:2));
    xspace = linspace(min_bb(1),max_bb(1),N);
    yspace = linspace(min_bb(2),max_bb(2),N);
    [X Y] = meshgrid(xspace,yspace);
    img1 = mk_image(i_gr.fwd_model,1);
    img1.calc_colours.npoints = N;
    M = calc_slices(img1,1);
    IN = M==1;
    n_px = nnz(IN);
    OUT = isnan(M);
%     Zpos=h*ones(1,n_px);
    %     xyzr = [X(IN)'; Y(IN)'; Zpos; r*ones(1,n_px); ];
    Z_layers=1;
    % order the X(IN) and Y(IN) according to their distances to point 0,0 
%     temp_X=0.8*X(IN);
%     temp_Y=0.8*Y(IN);
    temp_X=0.9*X(IN);
    temp_Y=0.9*Y(IN);
    [sort_,ind]=sort(temp_X.^2+temp_Y.^2);
    X_all=repmat(temp_X(ind),Z_layers,1);
    Y_all=repmat(temp_Y(ind),Z_layers,1);
    if Z_layers~=1
        hspace=linspace(h, max_bb(3),Z_layers);
        Z_all=hspace'*ones(1,n_px);
        Z_all=reshape(Z_all',size(Z_all,1)*size(Z_all,2),1);
    else 
        Z_all=h*ones(1,n_px)';
    end
    %     xyzr = [X_all'; Y_all'; Z_all'; r*ones(1,length(Z_all)); ];
    xyzr = [[0;X_all]'; [0;Y_all]'; [h;Z_all]'; r*ones(1,length(Z_all)+1); ]; %add a central point. assume x=0,y=0;
%     [xyzr]=del_out_map(img1,N,xyzr);
    Zpos=Z_all';
end

[vh_sim,vi_sim] = simulate_movement(img, xyzr);
imgr= inv_solve(i_gr,vh_sim,vi_sim);

% calculate the GREIT parameters
figure
levels =[inf,inf,Zpos];
imgr.show_slices.img_cols =10;
show_slices(imgr, levels);
imgr.calc_colours.npoints = 64;
imgr.calc_slices.levels=levels;
SimInLayer=length(xyzr);
params = eval_GREIT_fig_merit(imgr, xyzr,SimInLayer);
figure
p_names = {'AR','PE','RES','SD','RNG'};
for i=1:5; subplot(5,1,i);
    plot(params(i,:)); ylabel(p_names{i});
    mean_save(i)=mean(abs(params(i,:)));
    std_save(i)=std(abs(params(i,:)));
    
end

%% plot performance map / only suitable for test 3 and 4
% [vh,vi] = simulate_movement(img, xyzr);
% if Test == 6 % test 4 but sorted
%     IN5 = repmat(IN(ind),[1 1 5]);
% else
    IN5 = repmat(IN,[1 1 5]);
% end
IN5 = shiftdim(IN5,2);
p_img = nan([5 size(IN)]);
p_img(IN5) = params(:,1:538);  % for N=32
% p_img(IN5) = params(:,1:2245);  % for N=64
% p_img(IN5) = params(:,1:740);
% p_img(IN5) = params;
figure
imagesc(squeeze(p_img(1,:,:)));

