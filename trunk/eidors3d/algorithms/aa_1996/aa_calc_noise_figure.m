function hparam= aa_calc_noise_figure( inv_model );
% AA_CALC_NOISE_FIGURE
% hparam= aa_calc_noise_figure( inv_model );
% inv_model  => inverse model struct
%
% In order to use this function, it is necessary to specify
% inv_model.hyperparameter. has the following fields
% hpara.func         = 'aa_calc_noise_figure';
% hpara.noise_figure = VALUE
% hpara.tgt_elems    = vector of element numbers of contrast in centre
%
% The NF parameter is defined in Adler & Guardo (1996), as
%   measurements => z (Mx1), image elements => x (Nx1)
%   NF = SNR_z / SNR_x
% SNR_z = mean(z) / std(z) = mean(z) /sqrt( M trace(Rn) )
% SNR_x = mean(x) / std(x) = mean(x) /sqrt( N trance(ABRnB'A)
%   where Rn = sigma_n x inv(W) = the noise covariance, iW= inv(W)
%     and A  = diag matrix s.t. A_ii = area of element i
%
% NF = mean(z)/mean(ABz) * sqrt(N/M) * sqrt( trace(ABiWB'A) / trance(iW) )
%
% given a reconstructor x=R(z), then
%   ABz = A R(z)
%
% NOTE: SNR _should_ be defined in terms of power! This defines
%       it in terms of images amplitude
% NOTE: M and N are wrong here. They should be trace/M NOT trace*M

% $Id: aa_calc_noise_figure.m,v 1.3 2005-06-29 16:39:28 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

invv = calc_noise_figure( inv_model, 1e-4);
keyboard
% FIXME: this is a hack for now
hparam = 1e-4;

% simulate homg data and a small target in centre
function [h_data, c_data]= simulate_targets( fwd_model, ctr_elems)

   %Step 1: homogeneous image
   sigma= ones( size(fwd_model.elems,1) ,1);

   img= eidors_obj('image', 'homogeneous image', ...
                   'elem_data', sigma, ...
                   'fwd_model', fwd_model );
   h_data=fwd_solve( img );

   %Step 1: inhomogeneous image with contrast in centre
   delta = 1e-2;
   sigma(ctr_elems) = 1 + delta;
   img= eidors_obj('image', 'homogeneous image', ...
                   'elem_data', sigma, ...
                   'fwd_model', fwd_model );
   c_data=fwd_solve( img );

% calculate the noise figure for inv_model parameters
% based on the provided hyperparameter hp
function NF = calc_noise_figure( inv_model, hp)

   [h_data, c_data]= simulate_targets( inv_model.fwd_model, ...
        inv_model.hyperparameter.tgt_elems);

   % create a new inv_model, except with hyperparameter
   imdl_hp = rmfield(inv_model,'hyperparameter');
   imdl_hp.hyperparameter.value= hp;
   imdl_hp = eidors_obj( 'inv_model', imdl_hp );

   tgt_img= inv_solve( imdl_hp, h_data, c_data); 

   data_prior= calc_data_prior( inv_model );

   

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
