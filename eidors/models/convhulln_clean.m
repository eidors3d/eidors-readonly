function [K,V] = convhulln_clean(pts,p);
% CONVHULLN_CLEAN: run convhulln and catch errors

% (C) 2018 Bartlomiej Grychtol, Andy Adler
% License: GPL version 2 or 3
% $Id$

  if ischar(pts) && strcmp(pts,'UNIT_TEST'); do_unit_test; return; end


  K=0; V=0; dim = size(pts,2);

  if size(pts,1) < 3; return; end 
  % move points to origin (helps for small elements at
  % large coordinates
  ctr = mean(pts);
  pts = bsxfun(@minus,pts,ctr);
  scale = max(abs(pts(:)));

  if scale == 0; return; end  %when there's only one point

  % scale largest coordinate to 1 (helps with precision)
  pts = pts ./ scale;
  p.scale = scale;

  % force thorough search for initinal simplex and
  % supress precision warnings
  pts= uniquetol(pts,1e-13,'ByRows',true,'DataScale',1);
  % This won't catch all cases
  if size(pts,1)<=size(pts,2); return; end
  if any(std(pts)<1e-14); return; end
  [K,V] = call_convhulln(pts,p);

 % numerical issues may produce tiny negative volume
 % undo scaling
  V = scale^dim * max(V,0);

function [K,V] = call_convhulln(pts,p);
  K=0; V=0;
  dim = size(pts,2);

try
    [K, V] = convhulln(pts,{'Qt Pp Qs'});
catch
    %redo it with "Joggle", but set to zero if small
    [K, V] = convhulln(pts,{'Qt Pp Qs QJ'});
    if V<1e-8; V=0; end
end
   

function [K,V] = call_convhulln_old(pts,p);
  K=0; V=0;
  dim = size(pts,2);

try
       [K, V] = convhulln(pts,{'Qt Pp Qs'});
catch err
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
        u = uniquetol(pts*p.scale,1e-14,'ByRows',true,'DataScale', 1);
        ok = ok | (size(u,1) <= dim  );
        if ~ok; switch dim;
           case 2; ok = colinear_test2d(u);
           case 3; ok = colinear_test3d(pts*p.scale);
           otherwise; error('not 2D or 3D');
        end; end
  end
%    Save cases were errors called
%      load -mat CHP.mat ptsi;
%      ptsi{end+1} = pts;
%      save -mat CHP.mat ptsi; 
  if ~ok
     if      eidors_debug('query','mk_tet_c2f:convhulln');
        debug_plot_tet(p.fmdl,p.rmdl,p.tri_todo,p.t, p.pts)
        keyboard
     elseif  eidors_debug('query','mk_tri2tet_c2f:convhulln');
        debug_plot_tri2tet(p.fmdl,p.rmdl,p.v,p.t, p.bot, p.top, p.pts)
        keyboard
     else
        fprintf('\n');
        eidors_msg(['convhulln has thrown an error. (',err.message,')', ...
           'Enable "eidors_debug on convhulln_clean" and re-run to see a debug plot'],0);
        rethrow(err);
     end
  end
end

% test for colinearity in 2D
function ok = colinear_test2d(u,ok) 
   ok = false;
   cp = bsxfun(@minus, u(2:end,:), u(1,:));
   l = sqrt(sum(cp.^2,2));
   cp = bsxfun(@rdivide, cp, l);
   u = uniquetol(cp,1e-14,'ByRows',true,'DataScale',1);
   ok = ok | size(u,1) == 1; % co-linear points

% test for colinearity in 3D
function ok = colinear_test3d(pts);
   ok = false;
   u12 = uniquetol(pts(:,1:2),1e-14,'ByRows',true,'DataScale',1);
   cp = bsxfun(@minus, u12(2:end,1:2), u12(1,1:2));
   l = sqrt(sum(cp.^2,2));
   cp = bsxfun(@rdivide, cp, l);
   % counteract colinear vectors in different directions
   cp = abs(cp); 
   un = uniquetol(cp,1e-12,'ByRows',true,'DataScale',1);
   ok = ok | size(un,1) == 1; % co-linear points
   if ok; return; end

   % test if all points lie on the top or bottom caps
   top = max(pts(:,3));
   bot = min(pts(:,3));
   ok = ok | all(abs(pts(:,3) - top) < eps);
   ok = ok | all(abs(pts(:,3) - bot) < eps);

function debug_plot_tet(fmdl,rmdl,tri_todo,t, pts)
   clf
   tri.nodes = fmdl.nodes;
   vox.nodes = rmdl.nodes;
   tri.type = 'fwd_model';
   vox.type = 'fwd_model';
   vox.elems = rmdl.elems(v,:);
   vox.boundary = p.vox.elems;
   tri.elems = fmdl.elems(tri_todo(t),:);
   show_fem(vox)
   hold on
   h = show_fem(tri);
   set(h,'EdgeColor','b')
   pts = bsxfun(@plus,pts*scale,ctr);
   plot(pts(:,1),pts(:,2),'o');
   hold off
   axis auto

function debug_plot_tri2tet(fmdl,rmdl,v,t, bot, top, pts)
   clf
   tet.nodes = fmdl.nodes;
   tri.nodes = repmat(rmdl.nodes(rmdl.elems(v,:),:),2,1);
   tri.nodes(:,3) = [repmat(bot,3,1); repmat(top,3,1)];
   tri.elems = [ 1 2 5 4
                 2 3 6 5
                 3 1 4 6];
   tri.boundary = tri.elems;
   tet.type = 'fwd_model';
   tri.type = 'fwd_model';
   tet.elems = fmdl.elems(t,:);
   clf
   show_fem(tri)
   hold on
   h = show_fem(tet);
   set(h,'EdgeColor','b')
%    pts = bsxfun(@plus,pts*scale,ctr);
   plot3(pts(:,1),pts(:,2),pts(:,3),'o');
   hold off
   axis auto

function do_unit_test
   t1.type = 'fwd_model'; t1.elems = [1 2 3 4];
   t1.nodes = [0 0 0; 0 1 0; 1 0 0; 0 0 1]; t2 = t1;
   unit_test_cmp('A',mk_tet_c2f(t1,t2),    1,10*eps);
   t2.nodes(end,end) = -1;
   unit_test_cmp('B',mk_tet_c2f(t1,t2),    0,10*eps);

