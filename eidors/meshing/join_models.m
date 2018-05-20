function [fmdlo]= join_models(fmdl1, fmdl2, tol)
% JOIN_MODELS: Join two fmdl structures to create one
%
% [fmdlu]= join_models(fmdl1, fmdl2, tol)
% fmdlu is the union of fmdl1 and fmdl2 with common nodes removed.
% tol is the minimum distance between nodes to consider them unique
% tol (default 1e-6) if not provided. 
%
% join_models only works on the elems, nodes and electrode fields;
% any c2f, mat_idx, stimulation structures must be redefined
%

% (C) 2018 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fmdl1) && strcmp(fmdl1,'UNIT_TEST'); do_unit_test; return; end

if nargin==2; tol = 1e-5; end

fmdlo = do_join_models(fmdl1,fmdl2,tol);

function fmdlo = do_join_models(fmdl1,fmdl2,thresh);
   n1 = fmdl1.nodes;
   n2 = fmdl2.nodes;
   Dnodes = 0;
   nD = size(n1,2); %dimension
   for i=1:nD
      Dnodes = Dnodes + abs(bsxfun(@minus, n1(:,i), n2(:,i)'));
   end
   % idx of nodes in fmdl2 equal to nodes in fmdl1
   idx = find(any( Dnodes<thresh,1 ));

   nn1 = num_nodes(fmdl1);
   fmdlo = eidors_obj('fwd_model','joined model');
   fmdlo.nodes = [fmdl1.nodes;fmdl2.nodes];
   fmdlo.elems = [fmdl1.elems;fmdl2.elems+nn1];
   fmdlo.gnd_node = fmdl1.gnd_node;
   oidx = 1:num_nodes(fmdlo);
   nno = num_nodes(fmdlo);
   oidx(idx + nn1)= [];
   nidx = zeros(num_nodes(fmdlo),1);
   nidx(oidx) = 1:length(oidx);
   for i= idx(:)'
      ff = find(Dnodes(:,i)<thresh);
      if length(ff)~=1; error('Degenerate models with equal nodes'); end
      nidx(i+nn1) = ff;
   end
   fmdlo.nodes = fmdlo.nodes(oidx,:);
   fmdlo.elems = reshape(nidx(fmdlo.elems(:)),[],nD+1);
   if isfield(fmdl1,'electrode')
      fmdlo.electrode = fmdl1.electrode;
   else
      fmdlo.electrode = struct('nodes',{},'z_contact',{});
   end
   if isfield(fmdl2,'electrode')
      for i=1:length(fmdl2.electrode);
         eli = fmdl2.electrode(i);
         eli.nodes = nidx(eli.nodes + nn1);
         fmdlo.electrode(end+1) = eli;
      end
   end
   fmdlo.boundary = find_boundary(fmdlo.elems);



function do_unit_test
   subplot(221); [fmdl1,fmdl2]=do_unit_test_2D_mdls;

   fmdlo= join_models(fmdl1, fmdl2);
   do_testing('join_models-2D',fmdlo,fmdl1,fmdl2)
   subplot(222); show_fem(fmdlo,[0,1,1]); axis off




function do_testing(txt,fmdlo,fmdl1,fmdl2)
   unit_test_cmp([txt,'-01'], ...
       size(fmdl1.elems,1) + size(fmdl2.elems,1), ...
       size(fmdlo.elems,1) );
   unit_test_cmp([txt,'-02'], ...
       size(fmdl1.nodes,1) + size(fmdl2.nodes,1) >= ...
       size(fmdlo.nodes,1), true );

   if ~isfield(fmdl1,'electrode'); fmdl1.electrode = struct; end
   if ~isfield(fmdl2,'electrode'); fmdl2.electrode = struct; end
   unit_test_cmp([txt,'-02'], ...
       length(fmdl1.electrode) + length(fmdl2.electrode), ...
       length(fmdlo.electrode));

function [fmdl1,fmdl2]=do_unit_test_2D_mdls;
   imdl = mk_common_model('a2c0',8); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<0','x','y','z'));

   imdl = mk_common_model('a2c0',8); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<-0.4','x','y','z'));
   idx  = fmdl.nodes(:,1)<-0.25;
   fmdl.nodes(idx,1) = -0.25;
   fmdl.nodes(:,1) = fmdl.nodes(:,1) + 0.25;
   fmdl1 = crop_model(fmdl,inline('x> 1.0','x','y','z'));
   fmdl2 = fmdl;
   fmdl2.nodes(:,1) = -fmdl2.nodes(:,1);
   fmdl2 = crop_model(fmdl2,inline('x+y<-1.28','x','y','z'));
   fmdl2 = crop_model(fmdl2,inline('y<-0.95','x','y','z'));
   fmdl2 = crop_model(fmdl2,inline('x<-1.20','x','y','z'));
show_fem(fmdl1,[0,1,1]);
hold on; hh=show_fem(fmdl2,[0,1,1]); set(hh,'EdgeColor',[0,0,1]);
hold off; xlim([-1.3,1.3]); axis off
