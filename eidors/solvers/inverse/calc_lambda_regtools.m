function lambdas = calc_lambda_regtools(imdl, vh, vi, type, doPlot)
%% CALC_LAMBDA_REGTOOLS: Find optimal hyperparameter by the L-curve (LCC) 
% criterion or the generalized cross-validation (GCV).
%   lambdas = calc_lambda_regtools(imdl, vh, vi, type, doPlot);
%
% Output:
%   lambdas   - "optimal" hyperparameter(s) determined using LCC or GCV
%
% Input:
%   imdl      - inverse model (EIDORS struct)
%   vh        - homogenous voltage matrix (of size nVtg x 1)
%   vi        - inhomogenous voltage matrix (of size nVtg x nFrames) including noise(!)
%   type      - type of approach used, either:
%               'LCC' (default), the L-curve criterion
%               'GCV', generalized cross-validation
%   doPlot    - will enable plotting if set to true (default = false)
%
% Example:
%   calc_lambda_regtools('unit_test');  
%
% NOTE
%   if vi contains multiple frames the returned values will contain an
%   "optimal" hyperparameter for each frame. An appropriate lambda can then 
%   be determined from the average (e.g. median) of these values.
%
% See also: RTv4manual.pdf (please note that all page numbers listed
% correspond to the ones written in the upper right corner, the effective
% PDF page number will be += 2).
%
% Nomenclature: Jacobian J is A; Prior R (not RtR) is L; Voltage v is b
%
% Fabian Braun, December 2016
%
% CITATION_REQUEST:
% AUTHOR: P C Hansen
% TITLE: Regularization tools version 4.0 for Matlab 7.3.
% JOURNAL: Numerical algorithms
% YEAR: 2007
% VOL: 46
% NUM: 2
% PAGE: S189-194
% DOI: 10.1007/s11075-007-9136-9
%

% (C) 2016 Fabian Braun. License: GPL version 2 or version 3
% $Id$

citeme(mfilename);

%% unit testing?
if ischar(imdl) && strcmpi(imdl, 'unit_test')
   doUnitTest();
   return;
end


%% set default inputs
if ~exist('type', 'var') || isempty(type)
	type = 'LCC';
end
if ~exist('doPlot', 'var') || isempty(doPlot)
	doPlot = false;
end

%% check for existence of the regtools package
if exist('regudemo.m')==2  % file is already on path
% Do nothing. We're OK
elseif exist('./regtools', 'dir') %check if in current folder
   addpath('./regtools');
%%% What should this do?
elseif exist([fileparts(mfilename('fullpath')), filesep, 'regtools'])
   addpath([fileparts(mfilename('fullpath')), filesep, 'regtools'])
else
   error('Regtools are required but are not available, please download them from <a href="matlab: web http://www.mathworks.com/matlabcentral/fileexchange/52-regtools -browser">File Exchange</a> or <a href="matlab: web http://www2.compute.dtu.dk/~pcha/Regutools/ -browser">P.C. Hansen''s website</a> and store them in the subfolder called ''regtools''. In order to allow for a fast execution it is recommended to disable (uncomment) all plotting functions in l_cuve.m and gcv.m.');
end

% AA: 3feb2017: Please make changes so
% 1. we don't call get_RM
% 2. we call calc_R_prior
%     fix calc_R_prior so it does what you want
% 3. rename to calc_lambda_regtools
% 4. change tutorial to call new name
% 5. Make changes to mk_GREIT_model
% 6. Merge these changes into mainline (if it works)
%    OR: delete mainline and svn mv 

%% prepare imdl
imdlTmp = imdl;
imdlTmp.prior_use_fwd_not_rec = 0;  
% if isfield(imdl.fwd_model,'coarse2fine')
%     imdlTmp.fwd_model = rmfield(imdlTmp.fwd_model,'coarse2fine');
% end
% if isfield(imdl, 'rec_model') && isfield(imdl.rec_model,'coarse2fine')
%     imdlTmp.rec_model = rmfield(imdlTmp.rec_model,'coarse2fine');
% end
img_bkgnd = calc_jacobian_bkgnd(imdlTmp);
A = calc_jacobian(img_bkgnd);
W = calc_meas_icov(imdlTmp);
L = calc_R_prior(imdlTmp);   

LtL = calc_RtR_prior(imdlTmp);
LtL_ = L'*L;    
% assert(all(LtL_(:) - LtL(:) < 100*eps), 'Prior differs too much!');

% check that measurement covariance matrix W is identity
assert(isequal(W, speye(size(W))));


%% (IMPORTANT!) bring generalized to standart form (section 2.6 p.21 of RTv4manual.pdf)
% L-curves and of generalized and standard form are equal this is 
% because they have identical norms see (section 2.6.3 p.24 of RTv4manual.pdf)
[A_s, ~, ~] = std_form(A, L, nan(size(vh,1),1));  % as L is square b won't be affected, only A
% [A_s,b_s,L_p,K,M] = std_form(A,L,b);
% NOTE: We need it in standard form as l_curve and gcv routines only accept this
[U_s, s_s] = csvd(A_s);

%% Iterate through all frames to get a range of lambdas
nFrames = size(vi,2);
lambdas = nan(nFrames,1);

progress_msg('Calculating lambda for each frame:', 0, nFrames);

if doPlot
    figure(); 
end

data_width= max(num_frames(vi), num_frames(vh));
vi = filt_data( imdl, vi, data_width );
vh = filt_data( imdl, vh, data_width );
B = vh - vi;

for iFrame = 1:nFrames
    
    progress_msg(iFrame, nFrames);
    
    %% prepare differential data of current frame
    b = B(:,iFrame);
	
	switch(lower(type))
		case 'lcc'
			%% L-curve (see section 2.5 p.20 of RTv4manual.pdf)
			% calculate and plot continuous l-curve (documentation on p.83 of RTv4manual.pdf)
			lambdas(iFrame) = l_curve(U_s,s_s,b);

			% add my own l-curve to plot for validation purposes
			if doPlot && (iFrame == nFrames)
				lInit = imdl.hyperparameter.value;
				lams = flip(logspace(log10(lInit*1E-3), log10(lInit*1E3), 10));
				lams = [lams lInit];

				clear myrho myeta
				for i=1:length(lams);
					imdl.hyperparameter.value = lams(i);
					RM = get_RM( imdl );
					myrho(i) = (norm(A*(RM*b) - b));
					myeta(i) = (norm(L*(RM*b)));
				end

				% plot it
				hold on;
				loglog(myrho(1:end-1), myeta(1:end-1), 'ob');
				hold on;
				loglog(myrho(end), myeta(end), 'og');
			end
		case 'gcv'
			%% gcv (see p.37 of RTv4manual.pdf)
			% documentation on p.65 of RTv4manual.pdf
			lambdas(iFrame) = gcv(U_s,s_s,b);
			
			% plot my own GCV for validation purposes
			if doPlot && (iFrame == nFrames)
				lInit = imdl.hyperparameter.value;
				lams = flip(logspace(log10(lInit*1E-3), log10(lInit*1E3), 10));
				lams = [lams lInit];

				clear myG
				for i=1:length(lams);
					imdl.hyperparameter.value = lams(i);
					RM = get_RM( imdl );
					rho = (norm(A*(RM*b) - b))^2;
					myG(i) = rho / (trace(eye(size(RM,2)) - A*RM)^2);
				end

				% plot it
				hold on;
				loglog(lams(1:end-1), myG(1:end-1), 'ob');
				hold on;
				loglog(lams(end), myG(end), 'og');
			end
		otherwise
			error('type not supported!');
	end
end

progress_msg('Calculating lambda for each frame:', inf);

end

% TODO: this code really needs to be cleaned, but not before eidors 3.4
function nf= num_frames(d0)
   if isnumeric( d0 )
      nf= size(d0,2);
   elseif d0(1).type == 'data';
      nf= size( horzcat( d0(:).meas ), 2);
   else
      error('Problem calculating number of frames. Expecting numeric or data object');
   end
end

% test for existance of meas_select and filter data
function d2= filt_data(inv_model, d0, data_width )
   if ~isnumeric( d0 )
       % we probably have a 'data' object

       d1 = [];
       for i=1:length(d0)
          if strcmp( d0(i).type, 'data' )
              d1 = [d1, d0(i).meas];
          else
              error('expecting an object of type data');
          end
       end

   else
      % we have a matrix of data. Hope for the best
      d1 = d0;
   end

   d1= double(d1); % ensure we can do math on our object

   if isfield(inv_model.fwd_model,'meas_select') && ...
     ~isempty(inv_model.fwd_model.meas_select)
      % we have a meas_select parameter that isn []

      meas_select= inv_model.fwd_model.meas_select;
      if     size(d1,1) == length(meas_select)
         d2= d1(meas_select,:);
      elseif size(d1,1) == sum(meas_select==1)
         d2= d1;
      else
         error('inconsistent difference data: (%d ~= %d). Maybe check fwd_model.meas_select',  ...
               size(d1,1), length(meas_select));
      end
   else
      d2= d1;
   end

   if nargin==3 % expand to data width
      d2_width= size(d2,2);
      if d2_width == data_width
         % ok
      elseif d2_width == 1
         d2= d2(:,ones(1,data_width));
      else
         error('inconsistent difference data: (%d ~= %d)',  ...
               d2_width, data_width);
      end
   end
end

function doUnitTest()
% inspired by the tutorial mentioned below:
% http://eidors3d.sourceforge.net/tutorial/EIDORS_basics/tutorial110.shtml
%

% Load some data
load iirc_data_2006

stim = mk_stim_patterns(16,1,[0,1],[0,1],{'meas_current'},1);

for iRun = 1
% for iRun = [0 1]
    % Get a 2D image reconstruction model
    if iRun 
        % more advanced 3D model which includes coarse2fine mapping which
        % makes all crash
        fmdl = mk_library_model('adult_male_16el');
        fmdl.stimulation = stim;
        opts = [];
        opts.noise_figure = 0.5;
        imdl = mk_GN_model(fmdl, opts, []);
    else
        % simple one
        imdl= mk_common_model('c2c');
        imdl.fwd_model.stimulation = stim;
        imdl.fwd_model = rmfield( imdl.fwd_model, 'meas_select');
        imdl.RtR_prior = @prior_tikhonov;
    end

    % load the real data
    vi = real(v_rotate)/1e4; vh = real(v_reference)/1e4;
    % allow double precision, else we run into (unexplainable) problems
    vi = double(vi); vh = double(vh);   

    % get the hyperparameter value via L-curve
    figure
    lambdas_lcc = calc_lambda_regtools(imdl,vh,vi,'LCC',true);

    % get the hyperparameter value via GCV
    lambdas_gcv = calc_lambda_regtools(imdl,vh,vi,'GCV',true);

    % visualize
    FramesOfInterest = [10 35 60 85];
    fig = figure(1 + iRun);
    subplot(121);
    imdl.hyperparameter.value = median(lambdas_lcc);
    imgs_lcc = inv_solve(imdl, vh, vi(:,FramesOfInterest));
    imgs_lcc.show_slices.img_cols = 1;
    show_slices(imgs_lcc);
    title('L-curve');

    subplot(122);
    imdl.hyperparameter.value = median(lambdas_gcv);
    imgs_gcv = inv_solve(imdl, vh, vi(:,FramesOfInterest));
    imgs_gcv.show_slices.img_cols = 1;
    show_slices(imgs_gcv);
    title('GCV');

end

end
