function [C,A] = manifold_patches(F)
%% Taken from: http://www.alecjacobson.com/weblog/?p=3618
% Credits: Alec Jacobson
  % Compute connected components of facets connected by manifold edges.
  %
  % Inputs:
  %   F  #F by simplex-size list of facets
  % Outputs:
  %   C  #F list of component ids
  %

  % simplex size
  ss = size(F,2);
  assert(ss == 3);

  % List of all "half"-edges: 3*#F by 2
  allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  % Sort each row
  sortallE = sort(allE,2);
  % IC(i) tells us where to find sortallE(i,:) in uE: 
  % so that sortallE(i,:) = uE(IC(i),:)
  [uE,~,IC] = unique(sortallE,'rows');
  % uE2F(e,f) = 1 means face f is adjacent to unique edge e
  uE2F = sparse(IC(:),repmat(1:size(F,1),1,ss)',1);
  % kill non-manifold edges
  uE2F(sum(uE2F,2)>2,:) = 0;
  % Face-face Adjacency matrix
  A = uE2F'*uE2F;
  % All ones
  A = A>0;
  % Connected components are patches
  %C = components(A); % alternative to graphconncomp from matlab_bgl
  [~,C] = graphconncomp(A);
end