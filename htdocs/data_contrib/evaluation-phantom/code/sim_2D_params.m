function [DET, greit_para]= sim_2D_params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_2D_params calculates detectability values and GREIT measures
% for a single object moving from center to edge of the tank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

XPos = [0:0.2:0.8]; idx=0;
for ctrx= XPos; idx=idx+1;
    ctry = 0.0; brad = 0.1;
    
    extra={'ball',sprintf(['solid ball = cylinder(', ...
        '%f,%f,0;%f,%f,1;%f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.10;'], ...
        ctrx, ctry, ctrx, ctry, brad) };
    [fmdl, img] = create_2D_model(extra);
    ctr = interp_mesh(fmdl);
    ctr=(ctr(:,1)-ctrx).^2 + (ctr(:,2)-ctry).^2;
    
    % calculate detectability value
    DET(idx) = calc_detect_val(img, brad, ctr);
    
    % calculate simulated voltages and reconstruct images
    img.elem_data = 1 + 0.01*(ctr < brad^2);
    imgr(idx) = calc_volts_inv(img)
end
%
imgr(1).elem_data = [imgr(:).elem_data];
figure; show_slices(imgr(1));

%% calculate and plot GREIT parameters
figure
YPos=zeros(1,idx); ZPos=zeros(1,idx);
xyzr_pt = [XPos;YPos;ZPos;ones(1,idx)*brad];
greit_para = calc_greit_para(imgr, xyzr_pt, XPos);

%% plot the detectability parameters
figure
plot(XPos, DET/max(DET),'linewidth',2);
ylabel('Normalized Distinguishability {\it z}');
xlabel('Radius');
end

function [fmdl, img] = create_2D_model(extra)
% create 2D model

fmdl= ng_mk_cyl_models(0,[16],[0.1,0,0.02],extra);
img= mk_image(fmdl,1);
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
img.fwd_model.stimulation= stim;

end

function greit_para = calc_greit_para(imgr,xyzr_pt,XPos)
% calculate and plot the GREIT PARAMETERS
greit_para = eval_GREIT_fig_merit(imgr(1), xyzr_pt);
p_names = {'AR','PE','RES','SD','RNG'};
for i=1:5; subplot(5,1,i);
    plot(XPos,greit_para(i,:)); ylabel(p_names{i});
    if i<5; set(gca,'XTickLabel',[]);end
end
xlabel('Radius fraction');
end

function DET = calc_detect_val(img, brad, ctr)
% calculate the detectability parameters
vol= get_elem_volume(img.fwd_model);
contr = 1e-2;
indx = zeros(size(vol));
indx(ctr < brad^2)=+ contr;
J  = calc_jacobian( img );
Jr = J*indx;
DET = 1e6*sqrt( Jr'*Jr );
end

function imgr = calc_volts_inv(img)
% calculate simulated voltages and reconstruct images
imgh= img; imgh.elem_data(:) = 1;
% calculate voltages
vi = fwd_solve(img);
vh = fwd_solve(imgh);

imdl= mk_common_model('c2c2',16);
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
imdl.fwd_model.stimulation = stim;
imdl.fwd_model.meas_select = meas_sel;
imdl.hyperparameter.value = .001;
imgr = inv_solve(imdl,vh,vi);
end
