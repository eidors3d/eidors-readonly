% Mesh Correspondance $Id: square_mesh02.m,v 1.1 2008-03-28 16:01:28 aadler Exp $

% Create grid based on mesh points
nn= 16; nl= 1:nn+1;
h0pts= s_mdl.nodes([nl,(nl-1)*(nn+1)+1],:);
h1pts= s_mdl.nodes([nl + nn*(nn+1),nl*(nn+1)],:);

z_depth= .1*ones(2*(nn+1),1);;
% Add the third dimension
v00pts= [h0pts, -z_depth]; v01pts= [h0pts, +z_depth];
v10pts= [h1pts, -z_depth]; v11pts= [h1pts, +z_depth];

xpts= [v00pts(:,1),v01pts(:,1),v11pts(:,1),v10pts(:,1),v00pts(:,1)]';
ypts= [v00pts(:,2),v01pts(:,2),v11pts(:,2),v10pts(:,2),v00pts(:,2)]';
zpts= [v00pts(:,3),v01pts(:,3),v11pts(:,3),v10pts(:,3),v00pts(:,3)]';
subplot(121)
plot3(xpts,ypts,zpts,'b');

axis([-1.1,+1.1,-1.1,+1.1,-0.4,+0.4]);
view(-47,28); axis square
