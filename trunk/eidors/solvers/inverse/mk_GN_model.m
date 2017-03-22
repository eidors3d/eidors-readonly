function imdl = mk_GN_model(img, opt, lambda)
% MK_GN_MODEL: make EIDORS inverse models using the GREIT approach
%   imdl = mk_GN_model(img,  opt, lambda)
%
% Output: 
%   imdl      - GREIT inverse model
%
% Parameters:
%   img       - image of the forward model including stimulation pattern, etc.
%   options   - structure with fields:
%     imgsz         - [xsz ysz] reconstructed image size in pixels 
%                     (default: [32 32])
%     square_pixels - forces square pixels if 1 (default: 0)
%     RtRprior      - desired regularization prior:
%                     'laplace' (default), 'noser' or 'tikhonov'
%     noise_figure - the noise figure (NF) to achieve. 
%     meas_icov    - the inverse of the noise covariance matrix.
%   lambda    - set a fixed hyperparameter, if empty (default value) 
%               the hyperparameter is chosen acording to the desired noise figure
%
% See also MK_GREIT_MODEL (from which this script was adapted)
%

% (C) 2016 Fabian Braun. License: GPL version 2 or version 3
% $Id$

    %% do unit testing?
    if ischar(img) && strcmpi(img, 'unit_test')
        do_unit_test();
        return;
    end	

    %% parse and prepare inputs
    if ~exist('lambda', 'var')
      lambda = [];
    end
    if strcmp(img.type, 'fwd_model')
        img = mk_image(img, 1);
    end
    fmdl = img.fwd_model;
    imdl = select_imdl( fmdl,{'Basic GN dif'});
   
    %% parse and default options
    opt = parse_options(opt,fmdl,imdl);

    %% Calculate rec_model (if absent)
    if ~isfield(img,'rec_model');
        opt.do_coarse2fine = 0;  
        [imdl.rec_model, imdl.fwd_model] = mk_pixel_slice(fmdl, opt.target_plane, opt);
        imdl.rec_model.nodes(:,3) = []; % the third dimension complicated display
        % medical orientation: NO, DO flip the x-axis of the model beforehand 
        imdl.rec_model.mdl_slice_mapper.y_pts = fliplr(imdl.rec_model.mdl_slice_mapper.y_pts);
    else
        imdl.rec_model = img.rec_model;
    end


    %% create GaussNewton reconstruction matrix    
    imdl = select_imdl(imdl, {'Basic GN dif'});
    
    imdl.fwd_model.coarse2fine = mk_coarse_fine_mapping(imdl.fwd_model, imdl.rec_model);
    
    imdl.inv_solve.calc_solution_error = 0 ;
    
    if isfield(img, 'elem_data') 
      % non-homogenous
      imdl.jacobian_bkgnd.value = img.elem_data;        
    else
      imdl.jacobian_bkgnd.value = 1;
    end
    
    % assign the desired prior
    if strcmpi(opt.RtRprior, 'laplace');
        imdl.RtR_prior = @prior_laplace;  
    elseif strcmpi(opt.RtRprior, 'noser');
        imdl.RtR_prior = @prior_noser;  
		imdl.prior_use_fwd_not_rec = true;	% else we have issues with c2f mapping
    elseif strcmpi(opt.RtRprior, 'tikhonov');
        imdl.RtR_prior = @prior_tikhonov;  
    else
        error(['undefined prior: ', opt.RtRprior]);
    end 
    
    % assign inverse noise covariance matrix
    if ~isempty(opt.meas_icov)
        imdl.meas_icov = opt.meas_icov;
    end
    
    %% determine hyperparameter (either via noise figure or set manually)
    if ~isempty(opt.noise_figure)        
        %% take the four elements closest to the center (taken from select_imdl) 
        % BUT: take the elements at height of the BELT!
        xyz_elems = interp_mesh( imdl.fwd_model );
        ctr_elems = mean(xyz_elems, 1);
        ctr_elems(3) = opt.target_plane;  % at belt height!
        xyz_elems = xyz_elems - ones(size(xyz_elems,1),1)*ctr_elems;
        d_elems   = sqrt( sum( xyz_elems.^2, 2 ));
        [~, e_idx] = sort(d_elems);

        imdl.hyperparameter.tgt_elems = e_idx(1:4);
        imdl.hyperparameter.noise_figure = opt.noise_figure;

        sv_log = eidors_msg('log_level'); eidors_msg('log_level', 2);
        imdl.hyperparameter.value = choose_noise_figure( imdl );
        eidors_msg('log_level', sv_log);
    elseif ~isempty(opt.image_SNR)       
        imdl.hyperparameter.image_SNR = opt.image_SNR;
        imdl.hyperparameter.value = choose_image_SNR(imdl);
    else
        imdl.hyperparameter.value = lambda;
    end

end


function opt = parse_options(opt,fmdl,imdl)

    if ~isfield(opt, 'imgsz'),     opt.imgsz = [32 32]; end
    if ~isfield(opt, 'square_pixels')
        opt.square_pixels = 0;
    end
    % Allow imdl.rec_model to overwrite options.imgsz
    if isfield(imdl,'rec_model') && ~isempty(imdl.rec_model)
        % this assumes rec_model is a rectangular grid, as it should
        opt.imgsz(1) = numel(unique(imdl.rec_model.nodes(:,1)))-1;
        opt.imgsz(2) = numel(unique(imdl.rec_model.nodes(:,2)))-1;
    end
    
    if ~isfield(opt, 'noise_figure'), opt.noise_figure = []; end
    
    if ~isfield(opt, 'meas_icov'), opt.meas_icov = []; end
    
    if ~isfield(opt, 'image_SNR'), opt.image_SNR = []; end
    
    % prior type for regularization 
    if ~isfield(opt, 'RtRprior')
        opt.RtRprior = 'laplace'; 
    end
        
    % Calculate the position of the electrodes
    Nelecs = length(fmdl.electrode);
    for i=1:Nelecs
       enodesi = fmdl.electrode(i).nodes;
       elec_loc(i,:) = mean( fmdl.nodes( enodesi,:),1 );
    end
    opt.elec_loc = elec_loc;
    
    try
        opt.normalize = fmdl.normalize_measurements;
    catch 
        opt.normalize = 0;
        eidors_msg('mk_GN_model: fmdl.normalize_measurements not specified, assuming 0');
    end
    
    % electrode target planes
    if ~isfield(opt, 'target_plane')
          opt.target_plane = mean(elec_loc(:,3));
    else
        t = opt.target_plane;
        minnode = min(fmdl.nodes);
        maxnode = max(fmdl.nodes);
        if t<minnode(3) || t>maxnode(3)
            warning('options.target_plane is outside the model!');
            eidors_msg('mk_GN_model: Resorting to default target_plane');
            opt.target_plane = mean(elec_loc(:,3));
        end
    end
    
end


function do_unit_test()

    % Recosntruct with GREIT
    fmdl = mk_library_model('cylinder_16x1el_fine');
    fmdl.nodes = fmdl.nodes/15; % make radius 1;
    fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
    opt.noise_figure = 1.0;
    imdl_gr = mk_GREIT_model(fmdl, 0.2, [], opt);

    opt = struct();
    opt.noise_figure = 1.0; 
    imdl_gn_lap = mk_GN_model(fmdl, opt);

    opt = struct();
    opt.noise_figure = 1.0; 
    opt.RtRprior = 'tikhonov';
    imdl_gn_tik = mk_GN_model(fmdl, opt);

    test_performance( { imdl_gr, imdl_gn_lap, imdl_gn_tik}, fmdl );

end
