function gallery_fwd_model = mk_gallery(elec_posn,data_tomel,n_contours,factor,levels)
%% Constructs a FEM gallery geometry in 2D or 3D and returns the
%  corresponding EIDORS "fwd_model" structure
%
% elec_posn = 4-column matrix as for example:
%
%    EZG04_Ring1= [ ...
%         101   1.037202       0.000871      -0.465710  ... 
%         132   1.000000      -0.650000      -0.440000];    
%
% where:
%   col 1 = electrode number (arbitrary)
%   col 2 = x position of the electrode (m)
%   col 3 = y position of the electrode (m)
%   col 4 = z position of the electrode (m)
%
% Only columns 3 and 4 (i.e. y and z) are used in this function, and the x
%   direction is taken along the gallery axis.
%
% data_tomel = TOMEL data are given as a matrix with at least 8 colums as,
% for example:
%
%  Data_Ring1_Dec2004_Wen32_1= [ ...
%       1     223    101    104    102    103   100    0.260163 ...
%     160    2285    132    115    105    110   100    0.041110];
%
% where:
%   col 1 = measurement number
%   col 2 = time elapsed since begining of measurement sequence (in seconds)
%   col 3 = A electrode ( = positive current electrode)
%   col 4 = B electrode ( = negative current electrode)
%   col 5 = M electrode ( = positive potential electrode)
%   col 6 = N electrode ( = negative potential electrode)
%   col 7 = intensity (mA) of the current injected from A to B
%   col 8 = voltage measured V(M) - V(N) for the measurement pattern
%
% n_contours: number of contours to create for making the gallery
%   cross-section (e.g. n_contours = 8)
%
% factor: interpolation factor applied to the original contour described by
%   the electrode points. For instance, factor = 3 means that the interpolated
%   contour shall count 3 more points than the original one.
%
% levels = []               => constructs a 2D FEM geometry
%        = [x1, x2,..., xN] => constructs a 3D FEM geometry by extruding 
%                              the 2D FEM geometry on planes x = x1, x2,...
%
% gallery_fwd_model = matlab structure conforming to the fwd_model EIDORS
%   convention with the additional entries necessary in other functions:
%
%   params.misc.compute_CCandSS= 'y' => the CC and SS time-consuming
%      matrices will be computed the next time the forward-modeling function
%      "dg_fwd_solve_coarse" will be called
%   params.misc.n_elec= n_elec
%   params.misc.n_circles=n_contours
%   params.misc.factor=factor
%
%   Only for 3D FEM:
%   params.misc.map_2Dto3D= a matrix mapping the 2D elems to 3D elems
%   params.misc.n_elems_2D= the number of 2D elems forming the basic 2D
%      cross-section of the gallery.
%
% Method: FEM geometry is constructed through the following steps:
%   step 1 - interpolate (by factor) the initial inner contour of the gallery
%           formed by the electrode positions (yi,zi)
%   step 2 - create n_contours-1 dilated versions of the inner countour. The
%           dilation factor is such that the FEM mesh is almost square. This
%           may be controlled by changing the variable "alpha" in
%           mk_2D_gallery. Each contour counts the same number of points as the
%           inner countour.
%   step 3 - do the triangular meshing by connecting the points from nearby
%           countours. A total of n_contours-1 layers of trianges are so
%           obtained, each counting n_elec*factor*2 triangles
%   step 4 - only if levels is not empty, extrude the 2D FEM to obtain a 3D FEM
%           model.
%
%   mk_gallery create single point electrodes with no contact impedance
%
% Dominique Gibert, April 2007
%
%%
n_elec= size(elec_posn,1);
cross_section= [elec_posn(:,3),elec_posn(:,4)];
params= mk_geom_gallery(n_contours,levels,cross_section,factor);
params.stimulation= mk_stim_patterns_tomel(elec_posn,data_tomel);
if isempty(levels)
    params.name= 'EIDORS 2D FEM gallery geometry';
%     params.misc.n_elems_2D= size(params.elems,1);
else
    params.name= 'EIDORS 3D FEM gallery geometry';
    params.misc.map_2Dto3D= elems_2Dto3D(2*n_elec*factor*(n_contours-1),length(levels));
%     params.misc.n_elems_2D= size(params.misc.map_2Dto3D,1);
%     params.misc.n_elems_3D= size(params.elems,1);
end
params.solve=      'dg_fwd_solve';
params.system_mat= 'dg_calc_system_mat';
params.jacobian=   'dg_calc_jacobian';
params.normalize_measurements= 0;
params.misc.perm_sym= '{n}';
params.misc.compute_CCandSS= 'y';
params.misc.n_elec= n_elec;
params.misc.n_circles= n_contours;
params.misc.factor= factor;
gallery_fwd_model= eidors_obj('fwd_model', params);
end

function param= mk_geom_gallery(n_contours,levels,cross_section,factor)
%% Auxiliary Function to Construct a FEM gallery geometry in 2D or 3D
%
% n_contours: number of contours to create for making the gallery
%   cross-section (e.g. n_contours = 8)
%
% levels = []               => constructs a 2D FEM geometry
%        = [x1, x2,..., xN] => constructs a 3D FEM geometry by extruding 
%                              the 2D FEM geometry on planes x = x1, x2,...
%
% cross_section = 2-column matrix containing the (x,y) points of the electrode
% positions forming the inner cross_section of the gallery
%
% factor: interpolation factor applied to the original contour described by
%   the electrode points. For instance, factor = 3 means that the interpolated
%   contour shall count 3 more points than the original one.
%
% Dominique Gibert, April 2007
%
%%
n_elec= size(cross_section,1);
[elem2,node2,bdy,point_elec_nodes]= mk_2D_gallery(cross_section,n_contours,factor);
if isempty(levels) % 2D
    idx= (0:n_elec-1)*length(point_elec_nodes)/n_elec + 1;
    elec_nodes= point_elec_nodes(idx);
    param.nodes = node2';
    param.elems = elem2';
else  %3D
    [elem3,node3,bdy,point_elec_nodes] = mk_3D_gallery(elem2,node2,levels,bdy,point_elec_nodes);
    idx= (0:n_elec-1)*length(point_elec_nodes)/n_elec + 1;
    half_lev= ceil( length(levels)/2 );
    elec_nodes= point_elec_nodes( half_lev, idx );
    param.nodes = node3';
    param.elems = elem3';
    param.misc.model2d.nodes = node2';
    param.misc.model2d.elems = elem2';
end
param.name= sprintf('EIT FEM by mk_gallery with N=%d levs=%d',n_contours,length(levels));

param.boundary = bdy';
param.gnd_node = size(param.nodes,1); % ground node placed somewhere on outer boundary
param.electrode =  mk_electrodes(elec_nodes);
end

function [ELEM,NODE,bdy_nodes,point_elec_nodes] = mk_2D_gallery(cross_section,n_contours,factor)
%% Constructs a 2D FEM gallery geometry
%
% cross_section = 2-column matrix containing the (x,y) points of the electrode
% positions forming the inner cross_section of the gallery
%
% n_contours: number of contours to create for making the gallery
%   cross-section (e.g. n_contours = 8)
%
% factor: interpolation factor applied to the original contour described by
%   the electrode points. For instance, factor = 3 means that the interpolated
%   contour shall count 3 more points than the original one.
%
% Method: FEM geometry is constructed through the following steps:
%   step 1 - interpolate (by factor) the initial inner contour of the gallery
%           formed by the electrode positions
%   step 2 - create n_contours-1 dilated versions of the inner countour. The
%           dilation factor is such that the FEM mesh is almost square. This
%           may be controlled by changing the variable "alpha" in
%           mk_2D_gallery. Each contour counts the same number of points as the
%           inner countour.
%   step 3 - do the triangular meshing by connecting the points from nearby
%           countours. A total of n_contours-1 layers of trianges are so
%           obtained, each counting n_elec*factor*2 triangles
%
% Dominique Gibert, April 2007
%
%%
% Sets the origin of coordinates at the barycenter of the cross_section
x_crs= cross_section(:,1)-mean(cross_section(:,1));
y_crs= cross_section(:,2)-mean(cross_section(:,2));
n_crs= size(x_crs,1);
% Interpolate the x's and y's of intermediate points along cross_section
% First point is added at the end to close the cross-section
n_bdy= n_crs*factor;
x= interp1((0:n_crs)*factor,[x_crs(:)' x_crs(1)],(0:n_bdy-1));
y= interp1((0:n_crs)*factor,[y_crs(:)' y_crs(1)],(0:n_bdy-1));
[phi,rho]= cart2pol(x,y);
ELEM=[];
NODE= [];
% First: compute node positions by "expanding" the inner boundary
% Here alpha = 1 => cells are almost square
%      alpha > 1 => cells are elongated in the radial direction
alpha = 1.5;
for k=1:n_contours
    NODE= [NODE [rho.*cos(phi);rho.*sin(phi)]];
    rho= rho+alpha*mean(rho)*2*pi/n_bdy;
end
for k=1:n_contours-1
    % the n's are node numbers
    n1= k*n_bdy+1;
    n2= (k+1)*n_bdy;
    n3= k*n_bdy;
    n4= (k-1)*n_bdy+1;
    n5= n4+1;
    n6= n1+1;
    ELEM=[ELEM; n1,n2,n3; n1,n3,n4; n1,n4,n5; n1,n5,n6];
    for j=3:2:n_bdy
        n1= k*n_bdy+j;
        n2= n1-1;
        n6= n1+1;
        n4= (k-1)*n_bdy+j;
        n3= n4-1;
        n5= n4+1;
        ELEM=[ELEM; n1,n2,n3; n1,n3,n4; n1,n4,n5; n1,n5,n6];
    end
end
ELEM= ELEM';
bdy_nodes= [(1:n_bdy) (1:n_bdy)+(n_contours-1)*n_bdy];
bdy_nodes= [bdy_nodes;bdy_nodes(2:n_bdy),1,bdy_nodes(n_bdy+2:2*n_bdy),(n_contours-1)*n_bdy+1];
point_elec_nodes= (1:factor:n_bdy);
end

function [ELEM,NODE,BDY,elec_nodes] = mk_3D_gallery(elem0,node0,levels,bdy,elec_nodes0)
%% Auxiliary Function to Construct a 3D FEM gallery geometry by extruding a
%  2D FEM structure
%
% Example of typical usage (first build 2D FEM structure, then extrude to obtain
%                           the 3D FEM structure):
%
%   [elem,node,bdy,point_elec_nodes] = mk_2D_gallery(cross_section,n_contours,factor);
%   [elem,node,bdy,point_elec_nodes] = mk_3D_gallery(elem,node,levels,bdy,point_elec_nodes);
%
% Dominique Gibert, April 2007
%
%%
d= size(elem0,1);
n= size(node0,2);
e= size(elem0,2);
elem_odd= [elem0([1 1 2 3],:),elem0([1 3 2 3],:),elem0([3 2 1 2],:)];
elem_even=[elem0([2 1 2 3],:),elem0([2 3 1 3],:),elem0([3 2 1 1],:)];
NODE= [node0;levels(1)*ones(1,n)];
ELEM= [];
bdy1= [bdy;bdy(1,:)];
bdy2= [bdy;bdy(2,:)];
bl= size(bdy,2);
BDY = [];
ln= length(levels);
for k=2:ln
    NODE=[NODE  [node0; levels(k)*ones(1,n)] ];
    BDY= [BDY, ...
        bdy1 + [(k-1)*n*ones(2,bl); (k-2)*n*ones(1,bl)], ...
        bdy2 + [(k-2)*n*ones(2,bl); (k-1)*n*ones(1,bl)] ];
    if rem(k,2)==1
        elem= elem_odd;
    else
        elem= elem_even;
    end
    ELEM= [ELEM (elem + ...
        [[(k-1)*n*ones(1,e);(k-2)*n*ones(3,e)] ...
        [(k-1)*n*ones(2,e);(k-2)*n*ones(2,e)] ...
        [(k-1)*n*ones(3,e);(k-2)*n*ones(1,e)]] ) ];
end
% Add top and bottom boundary
BDY= [elem0,BDY,elem0+n*(ln-1)];
% elec_nodes is all nodes for all layers
elec_nodes= ones(ln,1)*elec_nodes0 + (0:ln-1)'*n*ones(1,length(elec_nodes0));
end

%param.electrode = mk_electrodes( elec_nodes );
% Create the electrode structure from elec_nodes
% Currently implements point electrodes with 
%   contact impedance of near zero.
function elec_struct = mk_electrodes(elec_nodes)
for i= 1:length( elec_nodes )
    elec_struct(i).nodes     = elec_nodes(i);
    elec_struct(i).z_contact = 0.001; % corresponds to 1 ohm
end
% Need to do it this way to be compatible accross versions
if ~exist('elec_struct');
    elec_struct= [];
end
end

