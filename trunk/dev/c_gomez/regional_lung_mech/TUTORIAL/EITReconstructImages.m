function [eitimages,s] = EITReconstructImages(eitdata)
%EITRECONSTRUCTIMAGES   Wrapper function to reconstruct images from the
%given EIT data and return a structured output containing the images and
%meta information.
%
% Signature:
% function [eitimages,s] = EITReconstructImages(eitdata)
%
% Input:
% eitdata   struct      scalar      .structure
%                                   .subject_ID
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .number_of_frames
%                                   .frame_rate
%                                   .voltage_data
%
% Output:
% eitimages struct      scalar      .structure
%                                   .subject_ID
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .image_series_length
%                                   .frame_rate
%                                   .difference_image_series
%                                   .absolute_iamge_series
%                                   .image_mask
%
% s         boolean     scalar      errors present: true
%
% Copyright C. Gomez-Laberge, November 2010.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 1
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITReconstructImages')
        display('Error EITReconstructImages: invalid arguments');
        error('Error EITReconstructImages: aborting execution');
end

% % Define function constants
ALGORITHM = 2; % 1=One step Gauss Newton with Laplace prior,
% 2=Electrode movement with Gaussian HPF prior.
K1 = 1;%0.95; % Shrink factor for reference voltage
K2 = 32; % pixels per image = K2^2

% Build FEM
[im,s] = BuildFem(eitdata.subject_type,eitdata.eit_device);

switch ALGORITHM
    case 1
        % Regularization hyperparameter chosen according to EIT system
        if strcmp(eitdata.eit_device,'Viasys')
            HP = 0.03;
        else % it's Draeger
            HP = 0.06;
        end
        im.hyperparameter.value = HP; % How shall we choose this?
        im.fwd_model.normalize_measurements = 1;
        im.RtR_prior = @laplace_image_prior;
    case 2
        % Regularization hyperparameter chosen according to EIT system
        if strcmp(eitdata.eit_device,'Viasys')
            HP = 0.06; P2= 1e-5;% Selected ad hoc minimize movement artefacts
%                 HP = 0.02; P2= 1e-5; % 1004 b2 c1
%                 HP = .2; P2= 1e-4; % 1006-RERV b2 c3 c4
%                 HP = 0.08; P2= 1e-4; % 1007-RERV c2
%                 HP = 0.2; P2= 1e-4; % 1007-RERV c3 d1
%                 HP = 0.10; P2= 1e-4; % 1008-RERV a b2 c1 c2
%                 HP = 0.10; P2= 1e-4; % 1008-RERV d5
           
        else % it's Draeger
            HP = 0.03; P2= 1e-5;% Selected ad hoc minimize movement artefacts
%                 HP = 0.16; P2= 1e-12;% 1011-RERV c1
%                 HP = 0.08; P2= 1e-12;% 1012-RERV

        end
        
        im.fwd_model.normalize_measurements=1;
        im.fwd_model.jacobian = @calc_move_jacobian;
        im.RtR_prior          = @elec_move_image_prior;
        %         im.elec_move_image_prior.RegC.func = @laplace_image_prior;
        im.elec_move_image_prior.RegC.func = @gaussian_HPF_prior;
        im.elec_move_image_prior.parameters = sqrt(P2/1);
        im.fwd_model.gaussian_HPF_prior.diam_frac = 0.2;
        im.hyperparameter.value = HP;
        im.inv_solve.select_parameters = 1:size(im.fwd_model.elems,1);
    otherwise
        display('Error FunctionTemplate: invalid reconstruction algorithm');
        s = true;
        error('Error FunctionTemplate: aborting execution');
end


% Calculated electrode resistances
%TO DO

% Define stimulation and measurement pattern
% if strcmp(eitdata.eit_device,'Draeger')
    [stim,meas,s] = BuildStimulationPattern(eitdata.eit_device);
    im.fwd_model.stimulation = stim;
    im.fwd_model.meas_select = meas;
% end

% Select difference image reconstruction algorithm
%TO DO of form: [imdiff,s] = CalculateDifferenceSolver(im)
imdiff = im;

% Reconstruct difference images
vd = eitdata.voltage_data;
if strcmp(eitdata.subject_type,'Test')
    [vref,s] = CalcReferenceImage(vd,'first10');
else
    [vref,s] = CalcReferenceImage(vd,'minim');
end
imgdelta = inv_solve(imdiff,vref*K1,vd);

% Set image size to K2 x K2 pixels
imgdelta.calc_colours.npoints = K2;
imgdelta.calc_slices.filter=ones(3)/9;
zdelta = -calc_slices(imgdelta);

% Calculate mask
mask = zeros(size(zdelta));
mask = ~isnan(zdelta(:,:,1));

% Scale zdelta change to match eitdata voltage change from reference
ztab = TabulateImageData(zdelta,mask);
znu = mean(ztab);
zmax = max(znu);
zmin = min(znu);
ztemp = zdelta-zmin;
nu = CalcGlobalSignal(eitdata,0);
scalefactor = 100*(max(nu)-min(nu))/min(nu);
zdelta = ztemp*scalefactor/zmax+zmin;

% Replace image NaN pixels with zeros
zdelta(isnan(zdelta)) = 0;


% % Select absolute image reconstruction algorithm
% %TO DO of form: [imabs,s] = CalculateAbsoluteSolver(im) % See EIDORS code below
% % imabs = im;
%
% Reconstruct absolute images
% % From EIDORS page:
% % In 3D, it's important to get the model diameter right, 2D is
% % imdl.fwd_model.nodes= imdl.fwd_model.nodes*15; % 30 cm diameter
% %
% % % Estimate the background conductivity
% % imgs = mk_image(imdl);
% % vs = fwd_solve(imgs); vs = vs.meas;
% %
% % pf = polyfit(vh,vs,1);
% %
% % imdl.jacobian_bkgnd.value = pf(1)*imdl.jacobian_bkgnd.value;
% % imdl.hyperparameter.value = imdl.hyperparameter.value/pf(1)^2;
% %
% % img = inv_solve(imdl, vh, vi);
z = [];

% Build eitimages structure
eitimages = struct('structure','eitimages',...
    'subject_id',eitdata.subject_id,...
    'subject_type',eitdata.subject_type,...
    'case_date',eitdata.case_date,...
    'file_name',eitdata.file_name,...
    'eit_device',eitdata.eit_device,...
    'maneuver',eitdata.maneuver,...
    'peep',eitdata.peep,...
    'deltap',eitdata.deltap,...
    'image_series_length',eitdata.number_of_frames,...
    'frame_rate',eitdata.frame_rate,...
    'difference_image_series',zdelta,...
    'absolute_image_series',z,...
    'image_mask',mask);

% End of function
end %function