function [K,V] = convhulln_clean(pts);
% CONVHULLN_CLEAN: run convhulln and catch errors

% (C) 2018 Bartlomiej Grychtol, Andy Adler
% License: GPL version 2 or 3
% $Id$

  if ischar(pts) && strcmp(pts,'UNIT_TEST'); do_unit_test; return; end

  K=0; V=0;

  dim = size(pts,2);

  if size(pts,1) < 3; return; end 
  % move points to origin (helps for small elements at
  % large coordinates
  ctr = mean(pts);
  pts = bsxfun(@minus,pts,ctr);
  scale = max(abs(pts(:)));

  if scale == 0; return; end  %when there's only one point

  % scale largest coordinate to 1 (helps with precision)
  pts = pts ./ scale;

  % force thorough search for initinal simplex and
  % supress precision warnings
  pts= uniquetol(pts,1e-13,'ByRows',true,'DataScale',1);
  % This won't catch all cases
  if size(pts,1)<=size(pts,2); return; end
  if any(std(pts)<1e-14); return; end
  [K,V] = call_convhulln(pts);

 % numerical issues may produce tiny negative volume
 % undo scaling
  V = scale^dim * max(V,0);
   

function [K,V] = call_convhulln(pts);
  DEBUG =  eidors_debug('query','convhulln_clean');
  dim = size(pts,2);

try
       [K, V] = convhulln(pts,{'Qt Pp Qs'});
catch err
%    Save cases were errors called
%      if ~exist('ok'); ptsi={}; end
%      ptsi{end+1} = pts;
%      save -mat CHP.mat ptsi; 
  ok = false;
  if exist('OCTAVE_VERSION')
     if strcmp(err.message,'convhulln: qhull failed')
        err.identifier =  'MATLAB:qhullmx:DegenerateData';
     end
        
  end
  switch err.identifier
     case {'MATLAB:qhullmx:DegenerateData', 'MATLAB:qhullmx:UndefinedError'}
        % border case point may be included multiple times.
        % this is OK... though we may miss cases where more
        % points should have been found but were not
        u = uniquetol(pts,1e-14,'ByRows',true,'DataScale', 1);
        ok = ok | (size(u,1) <= dim  );
        if ~ok; switch dim;
           case 2; ok = colin_test2d(u);
           case 3; ok = colin_test3d(u);
           otherwise; error('not 2D or 3D');
        end; end
  end
  if ~ok
     if DEBUG || eidors_debug('query','mk_tet_c2f:convhulln');
        tri.nodes = fmdl.nodes;
        vox.nodes = rmdl.nodes;
        tri.type = 'fwd_model';
        vox.type = 'fwd_model';
        vox.elems = rmdl.elems(v,:);
        vox.boundary = vox.elems;
        tri.elems = fmdl.elems(tri_todo(t),:);
        clf
        show_fem(vox)
        hold on
        h = show_fem(tri);
        set(h,'EdgeColor','b')
        pts = bsxfun(@plus,pts*scale,ctr);
        plot(pts(:,1),pts(:,2),'o');
        hold off
        axis auto
        keyboard
     else
        fprintf('\n');
        eidors_msg(['convhulln has thrown an error. (',err.message,')', ...
           'Enable eidors_debug on mk_tri_c2f and re-run to see a debug plot'],0);
        rethrow(err);
     end
  end
end

% test for colinearity in 2D
function ok = colin_test2d(u,ok) 
   ok = false;
   cp = bsxfun(@minus, u(2:end,:), u(1,:));
   l = sqrt(sum(cp.^2,2));
   cp = bsxfun(@rdivide, cp, l);
   u = uniquetol(cp,1e-14,'ByRows',true,'DataScale',1);
   ok = ok | size(u,1) == 1; % co-linear points

% test for colinearity in 3D
function ok = colin_test3d(u);
   ok = false;
   u12 = uniquetol(u(:,1:2),1e-14,'ByRows',true,'DataScale',1);
   cp = bsxfun(@minus, u12(2:end,1:2), u12(1,1:2));
   l = sqrt(sum(cp.^2,2));
   cp = bsxfun(@rdivide, cp, l);
   % counteract colinear vectors in different directions
   cp = abs(cp); 
   un = uniquetol(cp,1e-12,'ByRows',true,'DataScale',1);
   ok = ok | size(un,1) == 1; % co-linear points
   if ok; return; end

   % test if all points lie on the top or bottom caps
   ok = ok | all(abs(pts(:,3) - top) < eps);
   ok = ok | all(abs(pts(:,3) - bot) < eps);

function do_unit_test
   t1.type = 'fwd_model'; t1.elems = [1 2 3 4];
   t1.nodes = [0 0 0; 0 1 0; 1 0 0; 0 0 1]; t2 = t1;
   unit_test_cmp('A',mk_tet_c2f(t1,t2),    1,10*eps);
   t2.nodes(end,end) = -1;
   unit_test_cmp('B',mk_tet_c2f(t1,t2),    0,10*eps);

