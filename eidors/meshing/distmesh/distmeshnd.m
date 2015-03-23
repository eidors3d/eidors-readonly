function [p,t]=distmeshnd(fdist,fh,h,box,fix,varargin)
%DISTMESHND N-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
%
%      P:           Node positions (NxNDIM)
%      T:           Triangle indices (NTx(NDIM+1))
%      FDIST:       Distance function
%      FH:          Edge length function
%      H:           Smallest edge length
%      BOX:         Bounding box [xmin,xmax;ymin,ymax; ...] (NDIMx2)
%      FIX:         Fixed node positions (NFIXxNDIM)
%      FDISTPARAMS: Additional parameters passed to FDIST
%
%   Example: Unit ball
%      dim=3;
%      d=inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);
%
%   See also: DISTMESH2D, DELAUNAYN, TRIMESH, MESHDEMOND.

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
% This file is taken from the 3.2 Release by Per-Olof Persson,
% available under the GPL version 2 or any later version.   
% Source is from http://www-math.mit.edu/~persson/mesh/

global distmesh_do_graphics; % flag do decide if we do graphics
maxiter = 500;

dim=size(box,2);
ptol=.01;
ttol=.1;
L0mult=1+.4/2^(dim-1);
deltat=.10; % Stimulation Time step - Speed or slow simulation
geps=1e-1*h; deps=sqrt(eps)*h;

% 1. Create initial distribution in bounding box
if dim==1
  p=(box(1):h:box(2))';
else
  cbox=cell(1,dim);
  for ii=1:dim
    cbox{ii}=box(1,ii):h:box(2,ii);
  end
  pp=cell(1,dim);
  [pp{:}]=ndgrid(cbox{:});
  for ii=1:dim; pp{ii} = pp{ii} + h/2*(rand(size(pp{ii}))-0.5); end
  p=zeros(prod(size(pp{1})),dim);
  for ii=1:dim
    p(:,ii)=pp{ii}(:);
  end
end

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p,varargin{:})<geps,:);
r0=feval(fh,p,varargin{:});
p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
N=size(p,1);

count=0;
p0=inf;
while 1
  % 3. Retriangulation by Delaunay
  if max(sqrt(sum((p-p0).^2,2)))>ttol*h
    p0=p;
    t=delaunayn(p);
    pmid=zeros(size(t,1),dim);
    for ii=1:dim+1
      pmid=pmid+p(t(:,ii),:)/(dim+1);
    end
    t=t(feval(fdist,pmid,varargin{:})<geps,:); %BUGFIX -geps to geps
    % 4. Describe each edge by a unique pair of nodes
    pair=zeros(0,2);
    localpairs=nchoosek(1:dim+1,2);
    for ii=1:size(localpairs,1)
      pair=[pair;t(:,localpairs(ii,:))];
    end
    pair=unique(sort(pair,2),'rows');
    % 5. Graphical output of the current mesh
    if distmesh_do_graphics
       if dim==2
         hh=trimesh(t,p(:,1),p(:,2),zeros(N,1),'edgecolor','black');
         view(2),axis equal,axis off,drawnow
       elseif dim==3
         if mod(count,5)==0
           simpplot(p,t,'p(:,2)>0');
           title(['Retriangulation #',int2str(count)])
           drawnow
         end
      end
    else
      if count ==0
         opt.final_msg = '';
         progress_msg('distmeshng: retriangulation #:',0,0,opt);
      else
         progress_msg(count,0);
      end
    end
    count=count+1;
  end

  
  % 6. Move mesh points based on edge lengths L and forces F
  bars=p(pair(:,1),:)-p(pair(:,2),:);
  L=sqrt(sum(bars.^2,2));
  L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2,varargin{:});
  L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);
  F=max(L0-L,0);
  Fbar=[bars,-bars].*repmat(F./L,1,2*dim);
  dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...
                 ones(size(pair,1),1)*[1:dim,1:dim], ...
                 Fbar,N,dim));
  dp(1:size(fix,1),:)=0;
  p=p+deltat*dp;

  % 7. Bring outside points back to the boundary
  d=feval(fdist,p,varargin{:}); ix=d>0;
  gradd=zeros(sum(ix),dim);
  for ii=1:dim
    a=zeros(1,dim);
    a(ii)=deps;
    d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
    gradd(:,ii)=(d1x-d(ix))/deps;
  end
  p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;

  % 8. Termination criterion
  maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
  if maxdp<ptol*h; break; end
  % maximum iterations
  if count>maxiter; break; end
end
if ~distmesh_do_graphics
    progress_msg(Inf);
end

% final delaunayn just to make sure triangularization returned is good
t=delaunayn(p);
