function [content, list] = simplex_content(simp, nodes)
% SIMPLEX_CONTENT
%function [content, list] = simplex_content(simp, nodes)
%
% Find the 'content' of one or more simplex 'simp',
% given the simplex node coordinates 'node'.
%
% This function can return the volume contained in a
% set of tetrahedra, area covered by a set of triangles
% (or tetrahedra faces), or the length of a set of
% line segments (or triangle or tetrahedra edges).
% Collectively, volume, area, or length of a simplex is
% it's "content."
%
% The sum is returned in "content." The content of each
% simplex is returned in "list." Content is always positive.
%
% Care is taken here, to sum the smallest content first
% to ensure numeric stability.
%
% (C) 2015, Alistair Boyle

% See http://www.mathpages.com/home/kmath664/kmath664.htm
% See http://www.math.uci.edu/~chenlong/programming.html (iFEM)

if isstr(simp) && strcmp(simp, 'UNIT_TEST'); list = []; content = do_unit_test; return; end

nsimp = size(simp,2); % nodes per simplex
ndim = size(nodes,2); % number of dimensions for nodes
ndd = nchoosek(ndim, nsimp-1); % number of dimensions to sum over
% ... for the simple case (e.g. 3d volume, etc.) ndd = 1 (nsimp == ndim+1)
nddj = nchoosek(1:ndim,nsimp-1); % possible coordinate arrangements, if we're not doing it all at once

% we can calculate the content using the determinant
list = zeros(size(simp,1),1);
if ndd >= 1
%   V = 1 / 3! | 1 x1 y1 z1 |
%              | 1 x2 y2 z2 |
%              | 1 x3 y3 z3 |
%              | 1 x4 y4 z4 |
%   A = 1 / 2! | 1 x1 y1 |
%              | 1 x2 y2 |
%              | 1 x3 y3 |
%   L = 1 / 1! | 1 x1 |
%              | 1 x2 |
   % we start by scaling the nodes to be 0 --> +1
   if 0
     ns = size(simp,1); % simplex rows (how many simplices)
     un = unique(simp(:)); % unique nodes, in the simplices we care about
     nm  = min(nodes(un,:)); % nodal min coords
     node_scale = max(range(nodes(simp(:),:))); % node scaling
     nodes(un,:) = (nodes(un,:) - repmat(nm,length(un),1))/node_scale;
   else
     node_scale = 1;
   end

   %% OUTER LOOP: this is the 3D volume, 2D area, 1D length
   % now calculate the scaled element volumes
   scale = 1 / factorial(nsimp-1);
   for i = 1:size(simp,1)
      %% INNER LOOP: higher dimensional nodes than simplices... calculate the n-D hypotenuse
      detM=zeros(1,ndd);
      for j = 1:ndd
        M = [ ones(nsimp,1) nodes(simp(i,:),nddj(j,:)) ];
        detM(j)=det(M);
      end
      list(i) = sqrt(sum(sort(detM.^2)));
   end
else
   error('can''t yet handle degenerate simplices: nsimp = ndim+1... but we could: will you code me?');
end

content = sum(sort(abs(list))) * (scale * node_scale);
list = list * (scale * node_scale);

function ok = do_unit_test();
  ok = 1;
  % a simple line [0] to [1]
  lns = [ 1 2; 2 3]; lnn = [0:2]'/2;
  [lnl, lne]=simplex_content(lns,lnn);
  ok = match(ok, lnl, 1, 'line total length');
  ok = match(ok, lne, [1 1]'/2, 'line element lengths');
  % a simple line [0] to [1]
  lns = [ 1 2; 2 3]; lnn = [0:4]'*0.25;
  [lnl, lne]=simplex_content(lns,lnn);
  ok = match(ok, lnl, 0.5, 'line total length 0.5');
  ok = match(ok, lne, [1 1]'*0.25, 'line element lengths 0.5');
  lns = [ 1 2; 2 3]; lnn = [0:4]'*0.1; lnn(4) = 1;
  [lnl, lne]=simplex_content(lns,lnn);
  ok = match(ok, lnl, 0.2, 'line total length 0.2');
  ok = match(ok, lne, [1 1]'*0.1, 'line element lengths 0.2');
  % a simple box [0 0] to [1 1]
  boxs = [ 1 2 4; 2 3 4]; boxn = [0 0; 1 0; 1 1; 0 1];
  [boxa, boxe]=simplex_content(boxs,boxn);
  ok = match(ok, boxa, 1, 'box total area');
  ok = match(ok, boxe, [1 1]'/2, 'box element areas');
  % a simple cube [0 0 0] to [1 1 1]
  cubes = [1 2 3 7; 1 3 4 7; 1 5 6 7; 1 5 7 8; 1 2 7 6; 1 4 8 7];
  cuben = [0 0 0; 1 0 0; 1 1 0; 0 1 0; ...
           0 0 1; 1 0 1; 1 1 1; 0 1 1];
  [cubev, cubee]=simplex_content(cubes,cuben);
  ok = match(ok, cubev, 1, 'cube total volume');
  ok = match(ok, cubee, ones(6,1)/6, 'cube element volumes');
  %... now try stepping down a level: find_boundary, then go again
  % box boundary length
  boxf = find_boundary(boxs);
  [boxl, boxe]=simplex_content(boxf,boxn);
  ok = match(ok, boxl, 4, 'box total edge area');
  ok = match(ok, boxe, ones(4,1), 'cbox element length');
  % cube surface area
  cubef = find_boundary(cubes);
  [cubea, cubee]=simplex_content(cubef,cuben);
  ok = match(ok, cubea, 6, 'cube total surface area');
  ok = match(ok, cubee, ones(12,1)/2, 'cube element face area');

  if ok
    disp('simplex_content: PASS');
  else
    disp('simplex_content: FAIL');
  end

function ok= match( ok, pat1, pat2, descr)
  okl = (length(pat1(:)) == length(pat2(:))) && ...
        (all(pat1(:) == pat2(:)));
  if okl
    fprintf('simplex_content: pass - %s\n',descr);
  else
    fprintf('simplex_content: fail - %s\n',descr);
    if numel(pat1) < 15
      if (length(pat1(:)) == length(pat2(:)))
        disp([pat1(:) pat2(:)]);
      else
        pat1
        pat2
      end
    end
  end
  ok = okl && ok;
