function Reg= aa_calc_image_prior( inv_model );
% AA_CALC_IMAGE_PRIOR calculate image prior
% Reg= aa_calc_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   diam_frac= inv_model.image_prior.parameters(1) DEFAULT 0.1

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_calc_image_prior.m,v 1.5 2007-08-29 09:15:28 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );
if isfield(inv_model,'aa_calc_image_prior')
    diam_frac= inv_model.aa_calc_image_prior.parameters(1);
else
    diam_frac= 0.1;
end

Reg = calc_Gaussian_HPF( pp.NODE, pp.ELEM, diam_frac );

% Calculate Gaussian HP Filter as per Adler & Guardo 96
% parameter is diam_frac (normally 0.1)
function filt= calc_Gaussian_HPF( NODE, ELEM, diam_frac)
  e= size(ELEM, 2);
  np= 512; % number of interpolation points
% np= 128;
  taille=max(NODE')-min(NODE');

  xc= mean(reshape(NODE(1,ELEM(:)),3,e))'/taille(1);
  yc= mean(reshape(NODE(2,ELEM(:)),3,e))'/taille(2);
  
  [x y]=meshgrid( ...
      linspace( min(NODE(1,:)), max(NODE(1,:)) ,np ), ...
      linspace( min(NODE(2,:)), max(NODE(2,:)) ,np )  ); 
  v_yx= [-y(:) x(:)];
  o= ones(np*np,1);
  filt= zeros(e);
  tourne= [0 -1 1;1 0 -1;-1 1 0];

  for j= 1:e
    if ~rem(j,20);
       fprintf('.'); % if isoctave; fflush(stdout); end
    end
    xy= NODE(:,ELEM(:,j))';
    a= xy([2;3;1],1).*xy([3;1;2],2) ...
         -xy([3;1;2],1).*xy([2;3;1],2);
    aire=abs(sum(a));
    endr=find(y(:)<=max(xy(:,2)) & y(:)>=min(xy(:,2)) ...
            & x(:)<=max(xy(:,1)) & x(:)>=min(xy(:,1)) )';
    aa= sum(abs(ones(length(endr),1)*a' ...
         +v_yx(endr,:)*xy'*tourne)');
    endr( find( abs(1 - aa / aire) > 1e-8 ) )=[];
    ll=length(endr); endr=endr-1;

    % (rem(endr,np) corresponde a y
    % (floor(endr/np)) corresponde a x
    ym= ones(e,1)*(rem(endr,np)/(np-1) - .5) ...
        -yc*ones(1,ll);
    xm= ones(e,1)*(floor(endr/np)/(np-1) - .5) ...
        -xc*ones(1,ll);

    beta=2.769/diam_frac.^2;
%   filt(:,j)=-aire/2*beta/pi*mean(...
    filt(:,j)=-beta/pi*sum( exp(-beta*(ym.^2+xm.^2))')';
  end %for j=1:ELEM
% filt=filt/taille(1)/taille(2)+eye(e);
  filt=filt/np^2+eye(e);
  filt= sparse(filt.*(abs(filt)>.001)); 
