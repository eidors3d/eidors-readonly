function Reg= prior_gaussian_HPF( inv_model );
% PRIOR_GAUSSIAN_HPF calculate image prior
% Reg= prior_gaussian_HPF( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   diam_frac= inv_model.fwd_model.prior_gaussian_HPF.diam_frac DEFAULT 0.1
%
% CITATION_REQUEST:
% AUTHOR: A Adler & R Guardo
% YEAR: 1996
% TITLE: Electrical impedance tomography: regularized imaging and contrast
% detection 
% JOURNAL: IEEE transactions on medical imaging
% VOL: 15
% NUM: 2
% PAGE: 170–9
% LINK: http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=491418


% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

fwd_model= inv_model.fwd_model;
try 
    diam_frac= fwd_model.prior_gaussian_HPF.diam_frac;
catch
    diam_frac= 0.1;
end

copt.cache_obj= {fwd_model.nodes, fwd_model.elems, diam_frac};
copt.fstr = 'prior_gaussian_HPF';
if elem_dim(fwd_model) == 2;
   Reg = eidors_cache(@calc_Gaussian_HPF, {fwd_model, diam_frac}, copt );
else
   warning('prior_gaussian_HPF: not yet able to generate 3D models');
   Reg = prior_laplace( inv_model );
end


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

function do_unit_test
  imdl = mk_common_model('a2c0',16);
  RtR = prior_gaussian_HPF(imdl);
  tt=[0.562239752317943, -0.117068756722254, -0.025875127622824, -0.117068756722254;
     -0.117068756722254,  0.562239752317943, -0.117068756722254, -0.025875127622824;
     -0.025875127622824, -0.117068756722254,  0.562239752317943, -0.117068756722254;
     -0.117068756722254, -0.025875127622824, -0.117068756722254,  0.562239752317943];
  unit_test_cmp('a2c2 :1', RtR(1:4,1:4),tt,1e-10);

  imdl = mk_common_model('a3cr',16);
  RtR = prior_gaussian_HPF(imdl);  %NOTE: Fix required
  tt = [6    -2     0     0; -2     6     0     0;
        0     0     6    -2;  0     0    -2     6];
  unit_test_cmp('a3cr :1', RtR(1:4,1:4),tt,1e-10);

