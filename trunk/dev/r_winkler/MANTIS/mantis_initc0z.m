function [c0,z] = mantis_initc0z(VM,VC,stim,el,rimp,reco,verbose)
%MANTIS_INITC0Z initialize background conductivity c0 and contact impedances z
% Inputs:
% VM        - measured potentials
% VC        - simulated potentials for unit conductivity
% stim      - "stimulation" struct from EIDORS
% el        - struct containing electrode data:
%             el.P - number of electrode planes
%             el.Npp - number of electrodes per plane
%             el.Sizes - electrode surface areas
% rimp      - reference impedances used for computing VC (cell array)
% reco.zmin - lower conductivity limit
% reco.cmin - lower contact impedance limit
% reco.cmax - upper contact impedance limit
% verbose   - level of verbosity: 0 = 0ff,
%                               1 = timing,
%                               2 = important info + timing,
%                               3 = detailed info.
% Outputs:
% c0 - background estimate
% z  - contact impedances
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

if nargin < 5;                     rimp = [];                                      end
if nargin < 6 || isempty(reco);    reco.zmin=1e-5; reco.cmin=1e-4; reco.cmax=1e+4; end
if nargin < 7 || isempty(verbose); verbose = 0;                                    end

if ~isfield(reco,'zmin');  reco.zmin = 1e-5; end
if ~isfield(reco,'cmin');  reco.cmin = 1e-4; end
if ~isfield(reco,'cmax');  reco.cmax = 1e+4; end
if ~isfield(reco,'estZ');  reco.estZ = true; end
if reco.cmax <= reco.cmin; reco.cmax = 2*reco.cmin; end

%% for multiple planes: decompose to single planes, average individual estimates
if el.P > 1
  if verbose >= 3
    fprintf('Decompose %i planes for background and contact estimate.\n',el.P);
  end
  [stimr,VMr] = mantis_decplanes(stim,VM,el.Npp,el.P);
  [~,    VCr] = mantis_decplanes(stim,VC,el.Npp,el.P);

  c0 = zeros(el.P,1);
  z  = zeros(el.P*el.Npp,1);
  
  % recursive call of initial guesses
  for k=1:el.P
    kvec = (k-1)*el.Npp+1:k*el.Npp;
    elk.P = 1; elk.Npr = el.Npp; elk.Sizes = el.Sizes(kvec);
    
    if isempty(rimp)
      [c0(k),ztmp] = mantis_initc0z(VMr{k},VCr{k},stimr{k},elk,[],verbose);
    else
      [c0(k),ztmp] = mantis_initc0z(VMr{k},VCr{k},stimr{k},elk,rimp(kvec),reco,verbose);
    end
    z(kvec) = cell2mat(ztmp);
  end
  c0 = mean(c0);
  z = num2cell(z);
  return;
end

%% estimate for single planes
NOSERest = false; % do not use estimate without contact impedances...

if isempty(rimp) || ~reco.estZ; NOSERest = true; end

S = numel(stim);    % number of stimulations
N = zeros(S,1);     % number of electrodes
M = zeros(S,1);     % number of measurements
Mrank = zeros(S,1); % measurement rank

nind = 0;

for k=1:S
  M(k)     = size(stim(k).meas_pattern,1);
  N(k)     = size(stim(k).meas_pattern,2);
  tmpind   = nind+1:nind+M(k);
  nind     = nind+M(k);
  
  if any(abs(sum(stim(k).meas_pattern)) > 1e-10)
    stim(k).meas_pattern = stim(k).meas_pattern - ...
      (1/size(stim(k).meas_pattern,1))*repmat(sum(stim(k).meas_pattern),size(stim(k).meas_pattern,1),1);
    VM(tmpind) = VM(tmpind) - (1/length(tmpind))*sum(VM(tmpind));
    VC(tmpind) = VC(tmpind) - (1/length(tmpind))*sum(VC(tmpind));
  end
  
  Mrank(k) = rank(full(stim(k).meas_pattern));
end

if any(~(N(1)==N))
  error('Number of electrodes changes in measurements!');
else
  N=N(1);
end

if any(~(M(1)==M))                                                        % number of measurements per stimulation changes -> do NOSER estimate
  NOSERest = true;                                                        % (unlikely case. could be done better maybe...)
else
  M = M(1);
  if max(Mrank) < N-1 % all measurement ranks insufficient for contact impedance estimation
    NOSERest = true;
  else
    validS = find(Mrank >= N-1);
    if length(validS) < 2
      NOSERest = true;
    else
      % transform valid measurements to unit basis
      for k=1:S
        tmpind = (k-1)*M+1:k*M;
        if ismember(k,validS)
          if N > M; tolMM = 1e-8; else tolMM = min(svd(full(stim(k).meas_pattern)))+1e-8; end
          VM(tmpind) = pinv(full(stim(k).meas_pattern),tolMM)*VM(tmpind);
          VC(tmpind) = pinv(full(stim(k).meas_pattern),tolMM)*VC(tmpind);
        else
          VM(tmpind) = NaN; VC(tmpind) = NaN;
        end
      end
      S = length(validS);
      stim = stim(validS);
      VM = VM(~isnan(VM));
      VC = VC(~isnan(VC));      
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NOSERest || M~=N
  if verbose >= 2
    fprintf('Incomplete data. Using NOSER estimate and z=0.1*area.\n');
  end
  c0 = min(max(sum(VC(:).^2)/sum(VC(:).*VM(:)),reco.cmin),reco.cmax);
  z = num2cell(max(0.1*el.Sizes,reco.zmin));
  
else

  IM = zeros(N,S);
  for k=1:S; IM(:,k)=full(stim(k).stim_pattern); end                        % write currents as matrix

  VM = reshape(VM,N,S);
  VC = reshape(VC,N,S);
  
  A = sum((VC-repmat(rimp./(el.Sizes),1,S).*IM).*IM)';                    % compute "a" coeff. (see Winkler & Rieder 2015)
  B = sum(ones(N,S).*(IM.^2))';                                           % compute "b" coeff.
  C = sum(VM.*IM)';                                                       % compute "c" coeff.
  if cond([A B]) > 1e+10
    warning('Bad current pattern for finding c0 and z. Rotation invariant?');
    % try DCT transform to resolve rotation invariance
    I = dctmtx(N)';
    I = I(:,2:end);
    
    if N > S; tolIM = 1e-8; else tolIM = min(svd(IM))+1e-8; end

    VM = VM*pinv(IM,tolIM)*I;
    VC = VC*pinv(IM,tolIM)*I;
    IM = I;
    
    A = sum((VC-repmat(rimp./(el.Sizes),1,S).*IM).*IM)';                  % compute "a" coeff. (see Winkler & Rieder 2015)
    B = sum(ones(N,S).*(IM.^2))';                                         % compute "b" coeff.
    C = sum(VM.*IM)';                                                     % compute "c" coeff.
    
    if cond([A B]) > 1e+10
      warning('Could not resolve bad condition. Results may be inaccurate.');
    end
  end
  
  rz = [A B]\C;                                                           % least squares fit (see W & R 2015)
  
  if rz(2) < 0; rz(2) = 0; rz(1) = A\C; end                               % least squares fit failed (most likely contact impedance is close to 0). set to 0 and re-estimate background.
  
  z = num2cell(max(rz(2)*el.Sizes,reco.zmin));                            % avoid negative contacts
  c0 = min(max(1/rz(1),reco.cmin),reco.cmax);                             % avoid negative conductivity
  
end

end

