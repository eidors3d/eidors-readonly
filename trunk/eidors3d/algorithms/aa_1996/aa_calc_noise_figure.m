function hparam= aa_calc_noise_figure( inv_model );
% AA_CALC_NOISE_FIGURE
% hparam= aa_calc_noise_figure( inv_model );
% inv_model  => inverse model struct
%
% The NF parameter is defined in Adler & Guardo (1996), as
%   measurements => z (Mx1), image elements => x (Nx1)
%   NF = SNR_z / SNR_x
% SNR_z = mean(z) / var(z) = mean(z) /sqrt( M trace(Rn) )
% SNR_x = mean(x) / var(x) = mean(x) /sqrt( N trance(ABRnB'A)
%   where Rn = sigma_n x inv(W) = the noise covariance, iW= inv(W)
%     and A  = diag matrix s.t. A_ii = area of element i
%
% NF = mean(z)/mean(ABz) * sqrt(N/M) * sqrt( trace(ABiWB'A) / trance(iW) )
%
% NOTE: SNR _should_ be defined in terms of power! This defines
%       it in terms of images amplitude

% $Id: aa_calc_noise_figure.m,v 1.2 2005-06-27 15:55:46 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );

% FIXME: this is a hack for now
hparam = 1e-4;

return % ----------------------
sig= mean(DVV(:,1:4)')';
 n_var= 1 ./prob_dir( zeros(1,e) );
%n_var= zeros(size(sig));
H= sparse(1:m,1:m, n_var.^(-param(2)) );
  dd= DVV'*H*DVV;
  dx= DVV'*H;
  pond= filt'*filt;
if exist('Noise_Figs')==1 && ~Ignore_Noise_Figs
  bonH=find(Noise_Figs(:,3)==param(2) & ...
            Noise_Figs(:,4)==dpos );
else
  bonH=[];
  Noise_Figs=[];
end

if isempty(bonH)
disp('bonH==[]');
  lmax= 1e7*trace(dd)/trace(pond);
  Z=(dd+lmax*pond)\dx;
disp([sum(sig) sum(AIRE'*Z*sig) sqrt(n_var'*n_var) nfg]);
  nf= sum(sig)*sqrt(sum(sum( ((AIRE*n_var').*Z).^2 ))) / ...
       sum(AIRE'*Z*sig)/sqrt(n_var'*n_var)*nfg;
  Z=(dd+.1*lmax*pond)\dx;
  nf0= sum(sig)*sqrt(sum(sum( ((AIRE*n_var').*Z).^2 ))) / ...
       sum(AIRE'*Z*sig)/sqrt(n_var'*n_var)*nfg;
  if nf0<nf | nf<0
    disp('negative or non decreasing behaviour in NF');
    nf0=inf;
    lmax=1e-12*lmax;
    gap=10;
    while 1
      Z=(dd+lmax*pond)\dx;
      nf= sum(sig)*sqrt(sum(sum( ((AIRE*n_var').*Z).^2 ))) / ...
         sum(AIRE'*Z*sig)/sqrt(n_var'*n_var)*nfg;
      disp([log10(lmax) nf gap nf0])
      if nf>nf0 | nf<0
        lmax=lmax/gap;
        nf=nf0;
        if gap==10;
          gap= 1.25;
        else
          break;
        end        
      end
      nf0=nf;
      lmax=lmax*gap;
    end %while 1
    
  end

  fprintf('Choix=%s Filt=%d Minimum NF=%3.3f Lambda=%g\n', ...
             ChoiX, param(3)*100, nf, lmax );
  if isoctave; fflush(stdout); end
  Noise_Figs=[Noise_Figs; nf lmax param(2) dpos];
  bonH= size(Noise_Figs,1);
end %if any(Noise_Figs(:,3)==param(2))

if param(1)<min(Noise_Figs(bonH,1))
  disp('Noise Figure demande est trop petit');
  return
elseif any( abs(log(param(1)./Noise_Figs(bonH,1))) <.001 )
  ll= bonH(find( abs(log(param(1)./ ...
           Noise_Figs(bonH,1))) <.001 ));
  ll= Noise_Figs(ll(1),2);
  regul= ll*pond;
  Z=(dd+regul)\dx;
else
  nf=Noise_Figs(bonH,1:2);
  lmax= min([ nf( find( nf(:,1)<param(1) ), 2);1e10]);
  lmin= max([ nf( find( nf(:,1)>param(1) ), 2);1e-10]);
  ll= lmin;  nf= 1e16;

  iterat=0; 
  while abs(log(nf/param(1)))>.001
    iterat= iterat+1;
    if iterat>30
      disp('trop d''iterations');
      return;
    end

    if nf> param(1)
      lmin= ll;
    else
      lmax= ll;
    end; 
  
    ll=sqrt(lmax*lmin);
    Z=(dd+ll*pond)\dx;

    nf= sum(sig)*sqrt(sum(sum( ((AIRE*n_var').*Z).^2 ))) / ...
        sum(AIRE'*Z*sig)/sqrt(n_var'*n_var)*nfg;

    fprintf('Choix=%s Filt=%d NoiseFig=%3.3f Lambda=%g\n', ...
             ChoiX, param(3)*100, nf, ll );
    if isoctave; fflush(stdout); end
  end

  if ~Ignore_Noise_Figs
    Noise_Figs= [Noise_Figs; [nf ll param(2) dpos]];
    regul= ll*dd;
    if ~isoctave
       eval(['save ' fichier ' filt Noise_Figs'], ...
             'disp(''Error saving file'')'  );
    elseif ~exist(fichier)
       [filt_i, filt_j, filt_v, filt_nr, filt_nc]= spfind(filt);
       eval(['save -binary oct_',fichier, ...
             ' filt_i filt_j filt_v filt_nr filt_nc Noise_Figs'], ...
             'disp(''Error saving file'')'  );
    end
  end
end

% calculate nf with covariance
%   nf= sum(sig)*sqrt(sum( ((AIRE'*Z).*n_var').^2 )) / ...
%       sum(AIRE'*Z*sig)/sqrt(n_var'*n_var);
