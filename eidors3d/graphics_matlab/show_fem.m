function show_fem( mdl, sels )
% SHOW_FEM: show the EIDORS3D finite element model
% mdl is a EIDORS3D model structure
% $Id: show_fem.m,v 1.1 2004-07-16 14:10:14 aadler Exp $

set(gcf, 'Name', mdl.name);
trimesh(mdl.boundary, ...
        mdl.nodes(:,1), ...
        mdl.nodes(:,2), ...
        mdl.nodes(:,3) );
axis('image');
set(gcf,'Colormap',[0 0 0]);
hold on;

ee= mdl.boundary;
for e=1:length(mdl.electrode)
    % find elems on boundary attached to this electrode
    nn=ones(size(ee,1),1)*mdl.electrode(e).nodes;
    oo=ones(1,size(nn,2));
    ec=zeros(size(ee));
    for i=1:size(ec,2);
       ec(:,i) = any( (ee(:,i)*oo==nn)' )';
    end
    sels= find(all(ec'));

    for u=1:length(sels)
        paint_electrodes(sels(u),mdl.boundary, ...
                         mdl.nodes);
    end
end
hidden('off');

function paint_electrodes(sel,srf,vtx);
%function paint_electrodes(sel,srf,vtx);
%
%Auxilary function which plots the electrodes red at the boundaries.
%
%
%
% sel = The index of the electrode faces in the srf matrix
%       sel can be created by set_electrodes.m 
% srf = the boundary faces (triangles)
% vtx = The vertices matrix.


l = srf(sel,1); m = srf(sel,2); n = srf(sel,3);

Xs = [vtx(l,1);vtx(m,1);vtx(n,1)];
Ys = [vtx(l,2);vtx(m,2);vtx(n,2)];
Zs = [vtx(l,3);vtx(m,3);vtx(n,3)];

patch(Xs,Ys,Zs,'r');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003, A. Adler 2004
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
