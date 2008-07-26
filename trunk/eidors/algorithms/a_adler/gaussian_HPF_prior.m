function Reg= gaussian_HPF_prior( inv_model );
% GAUSSIAN_HPF_PRIOR calculate image prior
% Reg= gaussian_HPF_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   diam_frac= inv_model.fwd_model.gaussian_HPF_prior.diam_frac DEFAULT 0.1

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

fwd_model= inv_model.fwd_model;
try 
    diam_frac= fwd_model.gaussian_HPF_prior.diam_frac;
catch
    diam_frac= 0.1;
end

cache_test_obj= {fwd_model.nodes, fwd_model.elems, diam_frac};
Reg = eidors_obj('get-cache', cache_test_obj, 'gaussian_HPF_prior');
if ~isempty(Reg)
   eidors_msg('gaussian_HPF_prior: using cached value', 3);
   return
end

Reg = calc_Gaussian_HPF( fwd_model, diam_frac );

cache_test_obj= {fwd_model.nodes, fwd_model.elems};
eidors_obj('set-cache', cache_test_obj, 'gaussian_HPF_prior', Reg);
eidors_msg('gaussian_HPF_prior: setting cached value', 3);

% Calculate Gaussian HP Filter as per Adler & Guardo 96
% parameter is diam_frac (normally 0.1)
function filt= calc_Gaussian_HPF( fmdl, diam_frac)
  ELEM= fmdl.elems';
  NODE= fmdl.nodes';


  e= size(ELEM, 2);
  np= 128;
  [x,xc,y,yc] = interp_points(NODE,ELEM,np);

  v_yx= [-y,x];
  o= ones(np*np,1);
  filt= zeros(e);
  tourne= [0 -1 1;1 0 -1;-1 1 0];

  for j= 1:e
%   if ~rem(j,20); fprintf('.'); end
    xy= NODE(:,ELEM(:,j))';
    a= xy([2;3;1],1).*xy([3;1;2],2) ...
         -xy([3;1;2],1).*xy([2;3;1],2);
    aire=abs(sum(a));
    endr=find(y<=max(xy(:,2)) & y>=min(xy(:,2)) ...
            & x<=max(xy(:,1)) & x>=min(xy(:,1)) )';
    aa= sum(abs(ones(length(endr),1)*a' ...
         +v_yx(endr,:)*xy'*tourne)');
    endr( find( abs(1 - aa / aire) > 1e-8 ) )=[];
    ll=length(endr); endr=endr-1;

    yp= rem(endr,np)/(np-1) - .5; % (rem(endr,np) corresponde a y
    ym= ones(e,1)*yp -yc*ones(1,ll);
    xp= floor(endr/np)/(np-1) - .5; % (floor(endr/np)) corresponde a x
    xm= ones(e,1)*xp -xc*ones(1,ll);

    beta=2.769/diam_frac.^2;
%   filt(:,j)=-aire/2*beta/pi*mean(...
    filt(:,j)=-beta/pi*sum( exp(-beta*(ym.^2+xm.^2))')';
  end %for j=1:ELEM
% filt=filt/taille(1)/taille(2)+eye(e);
  filt=filt/np^2+eye(e);
  filt= ( filt+filt' )/ 2;
  filt= sparse(filt.*(abs(filt)>.003)); 

function [x,xc,y,yc] = interp_points(NODE,ELEM,np);
  taille=max(NODE')-min(NODE');
  e= size(ELEM, 2);

% Triangles of each shape
  xt= reshape(NODE(1,ELEM(:)),3,e)';
  yt= reshape(NODE(2,ELEM(:)),3,e)';

% We want center [1,1,1]/3 and edges [4,1,1]/6
  pts= [2,2,2;4,1,1;1,4,1;1,1,4]'/6;
  xp= xt*pts;
  yp= yt*pts;
  
  [x y]=meshgrid( ...
      linspace( min(NODE(1,:)), max(NODE(1,:)) ,np ), ...
      linspace( min(NODE(2,:)), max(NODE(2,:)) ,np )  ); 
% Add the basic interpolation points to those based on the
%  elements
  x= [x(:);xp(:)]; 
  y= [y(:);yp(:)]; 

  xc= mean(xt,2)/taille(1);
  yc= mean(yt,2)/taille(2);
