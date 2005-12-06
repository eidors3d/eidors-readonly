function [index_simp]=edge_refined_elem_mapper( mdl_coarse, mdl_dense)
% EDGE_REFINED_ELEM_MAPPER: map elements from coarse to dense model
% Calculates the index array mapping each dense mesh (from netgen)
%  simp onto a coarse mesh (from netgen) simp.
%
% Usage:
%  [index_simp]=edge_refined_elem_mapper( mdl_coarse, mdl_dense)
%
% (C) 2005 David Stephenson. Licensed under GPL v 2
% $Id: edge_refined_elem_mapper.m,v 1.3 2005-12-06 18:28:50 aadler Exp $

index_simp = eidors_obj('get-cache', mdl_dense, 'index_simp', mdl_coarse);
if ~isempty(index_simp)
    eidors_msg('edge_refined_elem_mapper: using cached value', 2);
    return
end

vtx_coarse =  mdl_coarse.nodes;
simp_coarse = mdl_coarse.elems;
vtx_dense =  mdl_dense.nodes;
simp_dense = mdl_dense.elems;

eidors_msg('edge_refined_elem_mapper: lookup_bld',4);
[lookup]=lookup_bld(simp_coarse,simp_dense);

eidors_msg('edge_refined_elem_mapper: center_of_simps',4);
[center_simp_dense]=center_of_simps(simp_dense, vtx_dense);

eidors_msg('edge_refined_elem_mapper: midpoints',4);
[vtx_midpoints_coarse]=midpoints(vtx_coarse,simp_coarse);

eidors_msg('edge_refined_elem_mapper: calc_h_refinement_centers',4);
[center_h_refined_simps]=calc_h_refinement_centers(simp_coarse,vtx_coarse,vtx_midpoints_coarse);

eidors_msg('edge_refined_elem_mapper: calc_simp_index',4);
[index_simp] = calc_simp_index(simp_dense,center_simp_dense,center_h_refined_simps,lookup);

% Cache the restult - it depends on both dense and coarse mdl
eidors_obj('set-cache', mdl_dense, 'index_simp', index_simp, mdl_coarse);
eidors_msg('inv_solve_dual_mesh: setting cached value', 2);


function [lookup]=lookup_bld(simp_coarse,simp_dense);

% This function calculates the lookup matrix of the dual mesh system (h-refinement) & is used for calculating index_simp. 
%
% simp_coarse = the coarse mesh simp array.
% simp_dense = the dense mesh simp array.

a=size(simp_coarse,1); % i.e the number of coarse mesh elements
b=size(simp_dense,1); % i.e the number of dense mesh elements

lookup=zeros(b,1);

x=1;
y=8;

for i=1:a;
    
    for ii=x:y;
        
        lookup(ii)=i;
        
        ii=ii+1;
        
    end
    
    x=x+8;
    y=y+8;
    
end


function [vtx_midpoints]=midpoints(vtx,simp);

% Calculates The x y z co-ordinates of the 6 midpoints of a 4 node simp.
%
% vtx = The vertices matrix
% simp = The simplices matrix
% vtx_midpoints = The x y z co-ordinates of the 6 midpoints of a 4 node simp

% Preallocation

vtx_midpoints=zeros(size(simp,1),18);

for i=1:size(simp,1);
    
    % Simp node 1~2
    x1=vtx(simp(i,1),1);
    y1=vtx(simp(i,1),2);
    z1=vtx(simp(i,1),3);
    x2=vtx(simp(i,2),1);
    y2=vtx(simp(i,2),2);
    z2=vtx(simp(i,2),3);

    
    x_mid=((x1+x2)/2);
    y_mid=((y1+y2)/2);
    z_mid=((z1+z2)/2);
    
    vtx_midpoints(i,1)=x_mid;
    vtx_midpoints(i,2)=y_mid;
    vtx_midpoints(i,3)=z_mid;
    
    % Simp node 1~3
    x1=vtx(simp(i,1),1);
    y1=vtx(simp(i,1),2);
    z1=vtx(simp(i,1),3);
    x2=vtx(simp(i,3),1);
    y2=vtx(simp(i,3),2);
    z2=vtx(simp(i,3),3);
    
    x_mid=((x1+x2)/2);
    y_mid=((y1+y2)/2);
    z_mid=((z1+z2)/2);
    
    vtx_midpoints(i,4)=x_mid;
    vtx_midpoints(i,5)=y_mid;
    vtx_midpoints(i,6)=z_mid;
    
    % Simp node 1~4
    x1=vtx(simp(i,1),1);
    y1=vtx(simp(i,1),2);
    z1=vtx(simp(i,1),3);
    x2=vtx(simp(i,4),1);
    y2=vtx(simp(i,4),2);
    z2=vtx(simp(i,4),3);
    
    x_mid=((x1+x2)/2);
    y_mid=((y1+y2)/2);
    z_mid=((z1+z2)/2);
    
    vtx_midpoints(i,7)=x_mid;
    vtx_midpoints(i,8)=y_mid;
    vtx_midpoints(i,9)=z_mid;
    
    % Simp node 2~3
    x1=vtx(simp(i,2),1);
    y1=vtx(simp(i,2),2);
    z1=vtx(simp(i,2),3);
    x2=vtx(simp(i,3),1);
    y2=vtx(simp(i,3),2);
    z2=vtx(simp(i,3),3);
   
    x_mid=((x1+x2)/2);
    y_mid=((y1+y2)/2);
    z_mid=((z1+z2)/2);
    
    vtx_midpoints(i,10)=x_mid;
    vtx_midpoints(i,11)=y_mid;
    vtx_midpoints(i,12)=z_mid;

    
    % Simp node 2~4
    x1=vtx(simp(i,2),1);
    y1=vtx(simp(i,2),2);
    z1=vtx(simp(i,2),3);
    x2=vtx(simp(i,4),1);
    y2=vtx(simp(i,4),2);
    z2=vtx(simp(i,4),3);
    
    x_mid=((x1+x2)/2);
    y_mid=((y1+y2)/2);
    z_mid=((z1+z2)/2);
    
    vtx_midpoints(i,13)=x_mid;
    vtx_midpoints(i,14)=y_mid;
    vtx_midpoints(i,15)=z_mid;
    
    % Simp node 3~4
    x1=vtx(simp(i,3),1);
    y1=vtx(simp(i,3),2);
    z1=vtx(simp(i,3),3);
    x2=vtx(simp(i,4),1);
    y2=vtx(simp(i,4),2);
    z2=vtx(simp(i,4),3);
    
    x_mid=((x1+x2)/2);
    y_mid=((y1+y2)/2);
    z_mid=((z1+z2)/2);
    
    vtx_midpoints(i,16)=x_mid;
    vtx_midpoints(i,17)=y_mid;
    vtx_midpoints(i,18)=z_mid;
    
    i=i+1;
    
end

function [center_h_refined_simps]=calc_h_refinement_centers(simp,vtx,vtx_midpoints);

% This function calculates the center of mass of each h-refined simp from the coarse mesh. 
%
% vtx      = The vertices matrix
% simp     = The simplices matrix
% vtx_midpoints  = The x y z co-ordinates of the 6 midpoints of a 4 node simp
% center_h_refined_simps = the center of mass of each h-refined simp from the coarse mesh

vtx_dave_1=[];
vtx_dave_2=[];
simp_dave=[1 5 6 7;2 5 8 9;3 6 8 10;4 7 9 10;5 6 7 9;5 6 8 9;6 7 9 10;6 8 9 10];
center_simp_proximity=[];

for i=1:size(simp,1);
    
	vtx_dave_1(1,1)=vtx(simp(i,1),1);
	vtx_dave_1(1,2)=vtx(simp(i,1),2);
	vtx_dave_1(1,3)=vtx(simp(i,1),3);
	
	vtx_dave_1(2,1)=vtx(simp(i,2),1);
	vtx_dave_1(2,2)=vtx(simp(i,2),2);
	vtx_dave_1(2,3)=vtx(simp(i,2),3);
	
	vtx_dave_1(3,1)=vtx(simp(i,3),1);
	vtx_dave_1(3,2)=vtx(simp(i,3),2);
	vtx_dave_1(3,3)=vtx(simp(i,3),3);
	
	vtx_dave_1(4,1)=vtx(simp(i,4),1);
	vtx_dave_1(4,2)=vtx(simp(i,4),2);
	vtx_dave_1(4,3)=vtx(simp(i,4),3);
	
	vtx_dave_2(1,1)=vtx_midpoints(i,1);
	vtx_dave_2(1,2)=vtx_midpoints(i,2);
	vtx_dave_2(1,3)=vtx_midpoints(i,3);
	
	vtx_dave_2(2,1)=vtx_midpoints(i,4);
	vtx_dave_2(2,2)=vtx_midpoints(i,5);
	vtx_dave_2(2,3)=vtx_midpoints(i,6);
	
	vtx_dave_2(3,1)=vtx_midpoints(i,7);
	vtx_dave_2(3,2)=vtx_midpoints(i,8);
	vtx_dave_2(3,3)=vtx_midpoints(i,9);
	
	vtx_dave_2(4,1)=vtx_midpoints(i,10);
	vtx_dave_2(4,2)=vtx_midpoints(i,11);
	vtx_dave_2(4,3)=vtx_midpoints(i,12);
	
	vtx_dave_2(5,1)=vtx_midpoints(i,13);
	vtx_dave_2(5,2)=vtx_midpoints(i,14);
	vtx_dave_2(5,3)=vtx_midpoints(i,15);
	
	vtx_dave_2(6,1)=vtx_midpoints(i,16);
	vtx_dave_2(6,2)=vtx_midpoints(i,17);
	vtx_dave_2(6,3)=vtx_midpoints(i,18);
	
	vtx_dave=[vtx_dave_1;vtx_dave_2];
	
	[center_simp_dave]=center_of_simps(simp_dave,vtx_dave);
	
	center_simp_proximity=[center_simp_proximity;center_simp_dave];
    
    vtx_dave_1=[];
    vtx_dave_2=[];
    
    i=i+1;
    
end

center_h_refined_simps=center_simp_proximity;



function [index_simp] = calc_simp_index(simp_dense,centre_simp_dense,center_h_refined_simps,lookup);

% Calculates the index array mapping each dense mesh (from netgen) simp onto a coarse mesh (from netgen) simp.
%
% simp_dense = The dense mesh simplices matrix
% centre_simp_dense = The center of mass of each (netgen) dense mesh element
% center_h_refined_simps = The center of mass of each h-refined simp from the coarse mesh
% lookup = the lookup matrix of the dual mesh system (h-refinement)

% Array pre-allocation

dist_simp=zeros(size(simp_dense,1),1);
index_simp=zeros(size(simp_dense,1),2);
mat_dense=zeros(size(simp_dense,1),1);

% Down to business ...

h = waitbar(0,'Calculating Simplex Map');

for id=1:size(simp_dense,1);   % for all dense center of simplicies

    waitbar(id/size(simp_dense,1))

    % find the x,y,z co-ord difference
    dx=centre_simp_dense(id,1)-center_h_refined_simps(:,1);
    dy=centre_simp_dense(id,2)-center_h_refined_simps(:,2);
    dz=centre_simp_dense(id,3)-center_h_refined_simps(:,3);
    
    dist_simp=sqrt((dx.^2)+(dy.^2)+(dz.^2));
       
    
    [m,I]=min(dist_simp);   % index out the minimum distance from the dense mesh to the id'th center of simplex
    
    index_simp(id,1)=lookup(I);

    index_simp(id,2)=m;   % write the actual minimum distance (as a quality control procedure)
    
end

close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) D.R Stephenson 2004
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version XXX
% MATLAB Version 6.5.0.180913a (R13)
% MATLAB License Number: 1560
% Operating System: Microsoft Windows XP Version 5.1 (Build 2600: Service Pack 1)
% Java VM Version: Java 1.3.1_01 with Sun Microsystems Inc. Java HotSpot(TM) Client VM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



