function [SNRmean, SE, debug] = calc_image_SNR(imdl, hyperparameter, doPlot)
%% CALC_IMAGE_SNR: Calculates the signal-to-noise ratio (SNR) in the image 
% domain as proposed by Braun et al. in:
% F Braun et al., A Versatile Noise Performance Metric for Electrical
% Impedance Tomography Algorithms, IEEE Trans. Biomed. Eng. 2017 (submitted).
%
%   [SNRmean, SE, debug] = calc_image_SNR(imdl, hyperparameter, doPlot)
%
% Output:
%   SNRmean   - mean of all SNR values as in equ. (9) in publication below
%   SE        - std of all SNR values as in equ. (9) in publication below
%   debug     - structure holding some information used for debug purposes
%
% Input:
%   imdl            - inverse model (EIDORS struct)
%      imdl.hyperparameter.roi_scaling_factor - amount of model shrinking
%                                   to determine ROI where seed n_T targets
%                                   (DEFAULT: 0.5);
%      imdl.hyperparameter.n_targets - number of targets to seed in ROI n_T
%                                   (DEFAULT: n_T = 200);
%      imdl.hyperparameter.target_radius - relative target radius r_T
%                                   (DEFAULT: r_T = 0.05);
%      imdl.hyperparameter.xyzr_targets - vector 4 x num targets n_T
%                                   specify targets manually [x y z r]
%   hyperparameter  - desired hyperparameter value, this will overwrite the
%                     imdl.hyperparameter.value and is for compatibility
%                     purposes with the function calc_noise_figure()
%   doPlot    will enable plotting if set to true (default = false)
%
%
% See also: CALC_NOISE_FIGURE
%
% Fabian Braun, December 2016
%
% CITATION_REQUEST:
% AUTHOR: F Braun et al.
% TITLE: A Versatile Noise Performance Metric for Electrical Impedance Tomography Algorithms
% JOURNAL: IEEE Transactions on Biomedical Engineering
% YEAR: 2017
% VOL: PP
% NUM: 99
% PAGE: 1
% DOI: 10.1109/TBME.2017.2659540
%

% (C) 2016 Fabian Braun. License: GPL version 2 or version 3
% $Id$

citeme(mfilename);

%% Configuration
% model shrunken by this factor delimits the ROI
ROI_SCALING_FACTOR = 0.5; 
% approximate number of targets to uniformly distribute in the ROI
N_DESIRED_TARGETS = 200;
% relative radius of the circular(2D)/spherical(3D) target (this is relative to the outer model radius) 
NORM_TGT_RADIUS = 0.05; % original
% TEZ region threshold 
TEZ_REGION_THRESH = 0.25;  % as ratio of maximum value : QUATER AMPLITUDE THRESHOLD


%% Parse and default inputs
if ~exist('doPlot', 'var') || isempty(doPlot)
    doPlot = false;
end
if isfield(imdl.hyperparameter, 'roi_scaling_factor') && ~isempty(imdl.hyperparameter.roi_scaling_factor)
    ROI_SCALING_FACTOR = imdl.hyperparameter.roi_scaling_factor;
end
assert(ROI_SCALING_FACTOR < 1, 'ROI must be scaled smaller than effective model');
if isfield(imdl.hyperparameter, 'n_targets') && ~isempty(imdl.hyperparameter.n_targets)
    N_DESIRED_TARGETS = imdl.hyperparameter.n_targets;
end
if isfield(imdl.hyperparameter, 'target_radius') && ~isempty(imdl.hyperparameter.target_radius)
    NORM_TGT_RADIUS = imdl.hyperparameter.target_radius;
end
if isfield(imdl.hyperparameter, 'xyzr_targets') && ~isempty(imdl.hyperparameter.xyzr_targets)
    xyzr_targets = imdl.hyperparameter.xyzr_targets;
else
    xyzr_targets = [];
end


%% input parsing
if nargin>=2 && numel(hyperparameter) == 1 && ~isempty(hyperparameter)
    imdl.hyperparameter.value = hyperparameter;
    % Remove function parameter because it will recurse
    try; imdl.hyperparameter = rmfield(imdl.hyperparameter,'func'); end
end


%% generate targets inside of rec_model
if isfield(imdl, 'rec_model')
    RecMdl = imdl.rec_model;   
else
    RecMdl = imdl.fwd_model;
end
MdlCtr = mean(RecMdl.nodes,1);

if isempty(xyzr_targets)
    % targets have not been defined, we'll generate them automatically 
    % by shrinking the model by the desired scaling factor and see the
    % descired number of targets inside the shrunk area
    %

    % first determine electrode level
    if isfield(RecMdl, 'mdl_slice_mapper')
        ElectrodeLevel = RecMdl.mdl_slice_mapper.level;    
    else
        if isfield(RecMdl.electrode(1), 'pos')
            elec_loc = cell2mat(cellfun(@(x) mean(x)', {RecMdl.electrode.pos}, 'uniformoutput', false))';
        elseif isfield(RecMdl.electrode(1), 'nodes')
            for i=1:length(RecMdl.electrode)
               enodesi = RecMdl.electrode(i).nodes;
               elec_loc(i,:) = mean( RecMdl.nodes( enodesi,:),1 );
            end  
        else
            error('not supported!');
        end
        ElectrodeLevel = mean(elec_loc,1);
    end
    
    % uniformly distribute N_DESIRED_TARGETS targets in the shrunken model 
    % at the level of the electrodes
    
    % first, shrink original model
    RecMdlShrunk = RecMdl;
    RecMdlShrunk.nodes = RecMdlShrunk.nodes - repmat(MdlCtr, size(RecMdlShrunk.nodes, 1), 1);
    RecMdlShrunk.nodes = RecMdlShrunk.nodes * ROI_SCALING_FACTOR;
    RecMdlShrunk.nodes = RecMdlShrunk.nodes + repmat(MdlCtr, size(RecMdlShrunk.nodes, 1), 1);

    % find boundary of shrunken model
    BoundaryNodes = find_boundary(RecMdlShrunk);
    % TODO: there's still a little problem here that we don't always
    % get a nicely oriented boundary, investigate this!
    Boundary = order_loop(RecMdlShrunk.nodes(BoundaryNodes, :));
    Boundary = [Boundary; Boundary(1,:)];
    % pp = fourier_fit(Boundary, 10); Boundary = fourier_fit(pp, linspace(0,1,20));

    % seed uniformly in rectangle and remove outliers
    AreaPoly = polyarea(Boundary(:,1), Boundary(:,2));
    Bounds = [max(Boundary); min(Boundary)];
    AreaRect = prod([Bounds(1,:) - Bounds(2,:)]);
    nUniformTgts = AreaRect * (N_DESIRED_TARGETS / AreaPoly);    

    % Size of the ROI in the two dimensions (x/y)...
    RoiSize = abs(diff(Bounds,[], 1));
    % ensure uniform spacing and scale according to differences in x/y size
    ScaleX = RoiSize(1) / RoiSize(2);
    ScaleY = 1./ScaleX;                 
    % create ROI centers
    [xx, yy] = meshgrid(linspace(Bounds(2,1), Bounds(1,1), ceil(ScaleX * sqrt(nUniformTgts))), ...
                        linspace(Bounds(2,2), Bounds(1,2), ceil(ScaleY * sqrt(nUniformTgts))));
    IsInside = inpolygon(xx(:), yy(:), Boundary(:,1), Boundary(:,2));
    nTgts = sum(IsInside);
    assert(abs((nTgts - N_DESIRED_TARGETS)/N_DESIRED_TARGETS) < 0.15, 'Cannot make desired number of targets');
    xx = xx(IsInside);
    yy = yy(IsInside);    
    zz = ones(size(xx))*ElectrodeLevel(3);
else
    % targets specified from outside the function, assign them properly
    xx = xyzr_targets(1,:)';
    yy = xyzr_targets(2,:)';
    zz = xyzr_targets(3,:)';
end

% set target size relative to maximal model radius
BoundsFull = [max(RecMdl.nodes); min(RecMdl.nodes)];
Rmodel = (max(BoundsFull(1,:) - BoundsFull(2,:)))/2;    % maximal model radius
Rtarget = Rmodel * NORM_TGT_RADIUS;
rr = ones(size(xx))*Rtarget;

if ~isempty(xyzr_targets)
    try
        rr = xyzr_targets(4,:)';    % overwrite target radii if existing
    end
end
    
%% generate differential voltages for each conductivity target 
img = mk_image(imdl.fwd_model, 1);
if elem_dim(imdl.fwd_model) == 3
    xyzr = [xx yy zz rr]';    
elseif elem_dim(imdl.fwd_model) == 2
    xyzr = [xx yy rr]';        
end
[vh, vi, xyzrOut, c2f] = simulate_movement(img, xyzr);
NotAssigned = ~ismember(xyzr', xyzrOut', 'rows');
assert(sum(NotAssigned) == 0, 'Error: target(s) got missing...');
vd = vi - repmat(vh, 1, size(vi,2));


%% get reconstruction matrix
RM = get_RM(imdl);
% calculate volume/area of each element in RecMdl 
RecMdlVols = get_elem_volume(RecMdl);


%% generate individual target evaluation zones (TEZs)
% first evaluate image response and take quater amplitude pixels as TEZ

% we calculate it the direct way as we can reuse it again further down!
imgrs = RM*vd;
imgrsNorm = imgrs;
imgrsNorm = imgrsNorm ./ repmat(max(imgrsNorm, [], 1), size(imgrsNorm,1), 1);
TEZs = double((imgrsNorm > TEZ_REGION_THRESH));

clear imgrsNorm;

ElemVols = spdiags(RecMdlVols, 0, length(RecMdlVols), length(RecMdlVols));

% correct ROI with element area/volume
TEZsNonNorm = ElemVols * TEZs; 
% normalize ROIs 
TEZs = TEZs ./ repmat(sum(TEZs,1), size(TEZs,1), 1);    


%% calculate target-wise distinguishability
% determine noise covariance matrix
if isfield(imdl.hyperparameter, 'SigmaN') && ~isempty(imdl.hyperparameter.SigmaN)
    SigmaN = imdl.hyperparameter.SigmaN;    
else
    SigmaN = speye(size(RM,2));		% assuming uniform and uncorrelated noise
end

Signal = diag(TEZsNonNorm' * imgrs);    % numerator in equation (9)
VarPerPixel = diag(RM*SigmaN*RM');      
Noise = sqrt(TEZs' * VarPerPixel);      % denominator in equation (9)    

% get the average of all SNRs
SNRs = Signal ./ Noise;                 
SNRmean = mean(SNRs); 
SE = std(SNRs);

if isnan(SNRmean)
    SNRmean = -inf;
elseif isinf(SNRmean)
    keyboard;
end

try
  eidors_msg('SNR = %e (hp=%e)', SNRmean, imdl.hyperparameter.value, 1);
end


%% debug info: fraction of amplitude response (AR) inside each TEZ
if nargout >= 3
    % target-wise amplitude response inside TEZ: in terms of energy
    ArInFrac = diag((TEZsNonNorm'*abs(RM*vd).^2) ./ (repmat(RecMdlVols, 1, size(TEZs,2))'*abs(RM*vd).^2));    
    debug.ArInFrac = ArInFrac;
    
    debug.SNRs = SNRs;
    debug.Signal = Signal;
    debug.Noise = Noise;
    debug.TEZs = TEZs;
    debug.TEZsNoNorm = TEZsNonNorm;
    
    debug.RoiBounds = Bounds;
    debug.RoiBoundary = Boundary;
    debug.MdlCtr = MdlCtr;
    debug.BoundsFull = BoundsFull;
    
    imgrsNorm = imgrs;
    imgrsNorm = imgrsNorm ./ repmat(max(imgrsNorm, [], 1), size(imgrsNorm,1), 1);
    debug.meanInNormImg = diag(TEZs'*imgrsNorm);
end


%% visualize for debug purposes
if doPlot
    fig = figure;
    set(fig, 'position', [  182         700        1531         313]);
    
    if isfield(imdl.fwd_model, 'coarse2fine')
        MapTgts2Img = imdl.fwd_model.coarse2fine'*c2f;
    else
        try
            imdl.fwd_model.coarse2fine = mk_coarse_fine_mapping(imdl.fwd_model, imdl.rec_model);
            MapTgts2Img = imdl.fwd_model.coarse2fine'*c2f;
        catch
            warning('Unable to make c2f mapping');
            MapTgts2Img = c2f;
        end
    end
    MapTgts2Img = MapTgts2Img ./ repmat(sum(MapTgts2Img,2), 1, size(MapTgts2Img, 2));
    MapTgts2Img(isnan(MapTgts2Img(:))) = 0;
    
    img = mk_image(RecMdl, nan);
    img.elem_data = MapTgts2Img * SNRs;
    SnrImg = calc_slices(img);
    
    img.elem_data = nan(size(img.elem_data));
    img.elem_data = MapTgts2Img * Signal;
    AmpSens = img.elem_data;
    AmpSensImg = calc_slices(img);
    
    img.elem_data = nan(size(img.elem_data));
    img.elem_data = MapTgts2Img * Noise;
    NoiseSens = img.elem_data;
    NoiseSensImg = calc_slices(img);
    
    sp1 = subplot(131);
    imagescnan(SnrImg);
    title(['SNR image: ', num2str(SNRmean, '%0.2d')]);
    colorbar; colormap jet;

    % plot signal sensitivity image
    sp2 = subplot(132);
    imagescnan(AmpSensImg);
    title(['Amplitude response: ', num2str(nanmean(AmpSens(:)), '%0.2d')]);
    colorbar;

    % plot noise sensitivity image
    sp3 = subplot(133);
    imagescnan(NoiseSensImg);
    title(['Noise sensitivity: ', num2str(nanmean(NoiseSens(:)), '%0.2d')]);
    colorbar;

    if doPlot < 2
        linkaxes([sp1 sp2 sp3], 'xy');
    end
   
    if doPlot == 2
       % some extra debugging 
       sp1 = subplot(131);
       imgTmp = img;
       iTgt = round(size(vd,2)/2);
       imgTmp.elem_data = RM*vd(:,iTgt);
       hh = show_fem(imgTmp);
       set(hh, 'edgecolor', 'none');
       hold on; 
       circle(xyzr(1:2, iTgt), xyzr(3, iTgt), 100, 'k');
       circle(xyzr(1:2, iTgt), Rtarget, 100, 'k');
       axis equal;
    end
    
    
    if doPlot == 3
       sp3 = subplot(133); 
       cla;
       hh = show_fem(RecMdl);
       set(hh, 'edgecolor', 'none')
       hold on;
       plot(Boundary(:,1), Boundary(:,2), '.-k');
       plot(xyzr(1,:), xyzr(2,:), '.r');
       plot(xyzr(1,:), xyzr(2,:), 'or');
    end
    
    if doPlot == 4
        % show pixel-wise noise sensitivity
       sp2 = subplot(132);
       cla;
       imgN = mk_image(imdl,0);
       imgN.elem_data = VarPerPixel;
       imgN.elem_data(sum(MapTgts2Img,2) == 0) = 0;
       imagescnan(calc_slices(imgN));
       title(['Pixel-wise NoiseSens: ', num2str(nanmean(imgN.elem_data(:)), '%0.2d')]);
       colorbar;        
    end
    
    if doPlot == 5
        sp1 = subplot(131);
        hist(SNRs); 
        xlim([0 max(xlim())]);
    end
    
end

end


