function [lambdas] = choose_lambda(imdl, vh, vi, type, doPlot)
%% Chooses the hyperparameter according to the L-curve (LCC) criterion 
% or the generalized cross-validation (GCV).
%
% Input:
%   - imdl      inverse model (EIDORS struct)
%   - vh        homogenous voltage matrix (of size nVtg x 1)
%   - vi        inhomogenous voltage matrix (of size nVtg x nFrames) including noise(!)
%   - type      type of approach used, either:
%               'LCC' (default), the L-curve criterion
%               'GCV', generalized cross-validation
%   - doPlot    will enable plotting if set to true (default = false)
%
% Output:
%   - lambdas   "optimal" hyperparameter(s) determined using LCC or GCV
%
% Note: if vi contains multiple frames the returned values will contain an
% "optimal" hyperparameter for each frame. An appropriate lambda can then 
% be determined from the average (e.g. median) of these values.
%
% See also: RTv4manual.pdf (please note that all page numbers listed
% correspond to the ones written in the upper right corner, the effective
% PDF page number will be += 2).
%
% Nomenclature: Jacobian J is A; Prior R (not RtR) is L; Voltage v is b
%
% Fabian Braun <fbn(ät)csem{dot}ch>, 17/12/2015
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
citeme(mfilename);

%% set default inputs
if ~exist('type', 'var') || isempty(type)
	type = 'LCC';
end
if ~exist('doPlot', 'var') || isempty(doPlot)
	doPlot = false;
end

%% check for existence of the regtools package
RegtoolsDir = [fileparts(mfilename('fullpath')), filesep, 'regtools'];
if exist(RegtoolsDir, 'dir')
	addpath(RegtoolsDir);
else
	error('Regtools are required but are not available, please download them from <a href="matlab: web http://www.mathworks.com/matlabcentral/fileexchange/52-regtools -browser">File Exchange</a> or <a href="matlab: web http://www2.compute.dtu.dk/~pcha/Regutools/ -browser">P.C. Hansen''s website</a> and store them in the subfolder called ''regtools''. In order to allow for a fast execution it is recommended to disable (uncomment) all plotting functions in l_cuve.m and gcv.m.');
end

%% prepare imdl
% A = calc_jacobian(img);
% LtL = calc_RtR_prior(imdl);
% W = calc_meas_icov(imdl);
[~, A, LtL, W] = get_RM( imdl );


% TODO: which method is best to get the R prior?
% (a) as defined in calc_R_prior (incomplete Cholesky factorization)
% L = ichol(LtL);     % we could do calc_R_prior but the coarse2fine bothers
% (b) doesn't give me what I want :-/
% L = calc_R_prior(imdl);    
% L=(L'*imdl.rec_model.coarse2fine)';
% (c) complete Cholesky fact. --> gives best results: L'*L closest to LtL
try
    L = chol(LtL);     % we could do calc_R_prior but the coarse2fine bothers
catch
    L = ichol(LtL);    % still unsure about that as L'*L does not lead to the same as LtL
end
LtL_ = L'*L;    
assert(all(LtL_(:) - LtL(:) < 100*eps), 'Prior differs too much!');

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

for iFrame = 1:nFrames
    
    progress_msg(iFrame, nFrames);
    
    %% prepare differential data of current frame
    b = vh - vi(:,iFrame);

	if doPlot
		figure(69); clf;
	end
	
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

