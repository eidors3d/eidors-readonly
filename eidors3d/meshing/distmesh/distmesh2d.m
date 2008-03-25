function [p,t]=distmesh2d(fd,fh,h0,bbox,pfix,varargin)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,PFIX,FPARAMS)
%
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      PFIX:      Fixed node positions (NFIXx2)
%      FPARAMS:   Additional parameters passed to FD and FH
%
%   Example: (Uniform Mesh on Unit Circle)
%      fd=inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%
%   Example: (Rectangle with circular hole, refined at circle boundary)
%      fd=inline('ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5))','p');
%      fh=inline('min(4*sqrt(sum(p.^2,2))-1,2)','p');
%      [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%
%   See also: MESHDEMO2D, DISTMESHND, DELAUNAYN, TRIMESH.

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
% This file is taken from the 3.2 Release by Per-Olof Persson,
% available under the GPL version 2 or any later version.   
% Source is from http://www-math.mit.edu/~persson/mesh/

dptol=.01; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0; deps=sqrt(eps)*h0;

% 1. Create initial distribution in bounding box (equilateral triangles)
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
p=[x(:),y(:)];                                       % List of node coordinates

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points
r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point
p=[pfix; p(rand(size(p,1),1)<r0./max(r0),:)];        % Rejection method

pold=inf;                                            % For first iteration
iteration=1;
while 1
  N=size(p,1);                                       % Number of points N
  % 3. Retriangulation by the Delaunay algorithm
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    t=delaunayn(p);                                  % List of triangles
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
    t=t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    [bars,i_bars,j_bars]=unique(bars,'rows');% Bars as node pairs
%   [bars,i_bars,j_bars]=unique(sort(bars,2),'rows');% Bars as node pairs
    % 5. Graphical output of the current mesh
    trimesh(t,p(:,1),p(:,2),zeros(N,1),'edgecolor','black')
    view(2),axis equal,axis off,drawnow
  end

  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths

  if 1
% Get length function at average
     hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
  else
% Get minimum length function at both ends
     hbars1=feval(fh,p(bars(:,1),:),varargin{:});
     hbars2=feval(fh,p(bars(:,2),:),varargin{:});
     hbars= min([hbars1,hbars2],[],2);
  end

  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths

  add_p= zeros(0,2);
if 0
  % 6.3 Desired length cannot be larger than tshape * connecting ones
  nt= size(t,1); % number of triangles
  min_L= zeros(nt,1);
  for i= 1:nt
     Lt = L(j_bars([i,i+nt,i+nt*2]));
     min_L(i) = min(Lt);
  
     if min(Lt)*5 < max(Lt); add_p= [add_p;pmid(i,:)];end
  end
  min_L = 2*[min_L;min_L;min_L];
  L0i= L0;
  L0 = min([L0,min_L(i_bars)],[],2);
end

  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)

  % 6.5. Points get sucked to the centroid to make them equilateral
if 0
  for i= 1:size(t,1)
     idx= t(i,:);
     dist_vtx = p(idx,:) - ones(3,1)*pmid(i,:);
     dist_mn  = sqrt(sum(dist_vtx.^2,2));
     dist_mn  = dist_mn/mean(dist_mn) - 1;
     dist_mn  = dist_mn.*(abs(dist_mn)>.4); %zero if too small
     F_adj    = (dist_mn*[1,1]) .* dist_vtx; 
     Fvec(idx,:) = Fvec(idx,:) + .3*F_adj;
  end
end

  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
  p=p+deltat*Ftot;                                   % Update node positions


  % 7. Bring outside points back to the boundary
  d=feval(fd,p,varargin{:}); ix=d>0;                 % Find points outside (d>0)
  dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
  dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; %    gradient
  p(ix,:)=p(ix,:)-[d(ix).*dgradx,d(ix).*dgrady];     % Project back to boundary

  % 8. Termination criterion: All interior nodes move less than dptol (scaled)
  if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end

  % 9. Add new points
  if size(add_p,1)>0 && rem(iteration,10)==1
     p= [p;add_p];
     pold=inf; 
  end
  iteration=iteration+1;
end
