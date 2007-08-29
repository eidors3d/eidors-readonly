function param= mk_circ_tank(rings, levels, elec_spec );
%MK_CIRC_TANK: make a cylindrical tank FEM geometry in 2D or 3D
% param= mk_circ_tank(rings, levels, elec_spec );
% 
% rings:  number of horizontal plane rings (divisible by 4)
% levels: vector of vertical placement of levels
%     for 2D mesh, levels = []
% 
% elec_spec: parameter to specify number of electrodes
%        specified as { 'opt1', val11, val12 , 'opt2', val21, val22 }
%
% elec_spec = scalar (divisible by 4)
%      - puts a single plane of electrodes in centre of cylinder
% eg. elec_spec  = 16
%
% elec_spec = { 'planes', n_elecs, elec_planes }
%      - puts plane each of n_elecs at planes specified by elec_planes
% eg. elec_spec  =  {'planes', 16, [2,6,8]}
%
% elec_spec = { 'zigzag', n_elecs, elec_planes }
%      - puts plane of n_elecs 'zigzagged' electrodes onto planes specified
%        1st elec on plane 2, 2nd elec on plane 6, 3rd on plane 2, etc 
% eg. elec_spec  =  {'zigzag', 16, [2,6]}
%      - Note, based on the restults of Graham et al (2006), zigzag
%        electrode placement is not recommended
%      - In order to implement the 'planar3d' pattern from this paper,
%        puts 2d electrodes onto rings ie [ ...  7  8  1  2  ...
%                                           ... 15 16  9 10  ... ]
%      ->use  elec_spec = { 'planes', n_elecs/2, elec_planes }
%
% mk_circ_tank creates simple, point electrodes. Improved models
%  may be created with create_tank_mesh_ng
%
% output:
%  param.name        Model name (if known) 
%  param.nodes       position of FEM nodes (Nodes x Dims) 
%  param.elems       definition of FEM elements (Elems x Dims+1) 
%  param.boundary    nodes of element faces on the medium surface 
%  param.gnd_node    Number of node connected to ground 
%  param.electrode   Vector (Num_elecs x 1) of electrode models (elec_model) 

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_circ_tank.m,v 1.19 2007-08-29 09:01:31 aadler Exp $

if rem(rings,4) ~= 0
   error('parameter rings and must be divisible by 4');
end

% parse easy case of electrode specifications
n_elec= [];
if size(elec_spec) == [1,1] if isnumeric(elec_spec)
   n_elec= elec_spec;
end; end

[elem, node, bdy, point_elec_nodes] = mk_2D_model( rings );

if isempty( levels ) % 2D
   
   if ~isempty( n_elec )
      idx= (0:n_elec-1)*length(point_elec_nodes)/n_elec + 1;
      elec_nodes= point_elec_nodes( idx );
   else
      error('2D models only support scalar electrode patterns');
   end
else  %3D
   [elem, node, bdy, point_elec_nodes] = mk_3D_model( elem, node, ...
                  levels, bdy, point_elec_nodes );

   if ~isempty( n_elec )
      idx= (0:n_elec-1)*length(point_elec_nodes)/n_elec + 1;
      half_lev = ceil( length(levels)/2 );
      elec_nodes= point_elec_nodes( half_lev, idx );
   else
      elec_nodes= electrode_pattern( point_elec_nodes, elec_spec );
   end

end

param.name= sprintf('EIT FEM by mk_circ_tank with N=%d levs=%d', ...
                    rings, length(levels) );
param.nodes = node';
param.elems = elem';
param.boundary = bdy';
param.gnd_node = 1; % node at bottom and center of the tank
param.electrode =  mk_electrodes( elec_nodes );

return;

% parse the elec_spec parameter
% elec_spec = { 'planes', n_elecs, elec_planes }
%      - puts plane each of n_elecs at planes specified by elec_planes
% eg. elec_spec  =  {'planes', 16, [2,6,8]}
%
% elec_spec = { 'zigzag', n_elecs, elec_planes }
%      - puts plane of n_elecs 'zigzagged' electrodes onto planes specified
%        1st elec on plane 2, 2nd elec on plane 6, 3rd on plane 2, etc 
% eg. elec_spec  =  {'zigzag', 16, [2,6]}
function elec_nodes= electrode_pattern( point_elec_nodes, elec_spec )
   elec_nodes= [];
   lpe = size(point_elec_nodes,2);
   nlev= size(point_elec_nodes,1);
   for i=1:3:length(elec_spec)-2
      spec = elec_spec{i};
      if      strcmp( spec,'planes' )
          n_elec= elec_spec{i+1};
          levs =  elec_spec{i+2};

          eidx= (0:n_elec-1);
          idx= round(eidx*lpe/n_elec) + 1;
          nodes= point_elec_nodes( levs, idx )';
          elec_nodes= [ elec_nodes; nodes(:) ];
      elseif  strcmp( spec,'zigzag' )
          n_elec= elec_spec{i+1};
          levs =  elec_spec{i+2};
          if any(levs > size(point_elec_nodes,1))
             error('requested electrode plane larger than FEModel');
          end

          eidx= (0:n_elec-1);
          idx= round(eidx*lpe/n_elec)*nlev + ...
               levs( rem( eidx, length(levs))+1);
          nodes= point_elec_nodes( idx );
          elec_nodes= [ elec_nodes; nodes(:) ];
      else
        error('elec_spec parameter not understood');
      end
   end

% Create a simple 2D regular mesh, based on N circular rings
%   and n_elec electrodes
function [ELEM, NODE, bdy_nodes, point_elec_nodes] = mk_2D_model( N );
  ELEM=[];
  NODE= [0;0];
  int=1;
  for k=1:N
    phi= (0:4*k-1)*pi/2/k;
    NODE= [NODE k/N*[sin(phi);cos(phi)]];

    ext= 2*(k*k-k+1);
    idxe=[0:k-1; 1:k];
    idxi=[0:k-1]; 
    elem= [ ext+idxe, ext+2*k+[-idxe idxe], ...
                     ext+rem(4*k-idxe,4*k), ...
            ext+idxe, ext+2*k+[-idxe idxe], ...
                     ext+rem(4*k-idxe,4*k);
            int+idxi, int+2*(k-1)+[-idxi, idxi], ... 
            int+rem(4*(k-1)-idxi, 4*(k-1)+(k==1) ) ...
            ext+4*k+1+idxi, ...
            ext+6*k+ [1-idxi 3+idxi], ...
            ext+8*k+3-idxi ];
    for j=1:k
      r1= rem(j+k-1,3)+1;
      r2= rem(j+k,3)+1;
      r3= 6-r1-r2;
      elem([r1 r2 r3],j+k*(0:7) )= elem(:,j+k*(0:7));
    end

    ELEM=[ ELEM elem(:,1:(8-4*(k==N))*k) ];
    int=ext;
  end %for k=1:N

  bdy_nodes= [ (ext  :ext+N*4-1) ; ...
               (ext+1:ext+N*4-1), ext ];
  point_elec_nodes= (ext):(ext+N*4-1) ;
 

% 'extrude' a 2D model defined by ELEM and NODE into a 3D model
% levels are defined by 'niveaux', 
% 2D parameters are ELEM, NODE, and bdy
function [ELEM, NODE, BDY, elec_nodes] = mk_3D_model( ...
     elem0, node0, niveaux, bdy, elec_nodes0 );

  d= size(elem0,1);       %dimentions+1
  n= size(node0,2);       %NODEs
  e= size(elem0,2);       %ELEMents     

  elem_odd= [elem0([1 1 2 3],:), ...
             elem0([1 3 2 3],:), ...
             elem0([3 2 1 2],:)]; 
  elem_even=[elem0([2 1 2 3],:), ...
             elem0([2 3 1 3],:), ...
             elem0([3 2 1 1],:)]; 
         
         
  NODE= [node0; niveaux(1)*ones(1,n) ];
  ELEM= [];
  bdy1= [bdy;bdy(1,:)];
  bdy2= [bdy;bdy(2,:)];
  bl= size(bdy,2);
  BDY = [];
 
  ln= length(niveaux);
  for k=2:ln
    NODE=[NODE  [node0; niveaux(k)*ones(1,n)] ];
    if rem(k,2)==1
        elem= elem_odd;
    else
        elem= elem_even;
    end
    ELEM= [ELEM (elem + ...
       [[(k-1)*n*ones(1,e);(k-2)*n*ones(3,e)] ...
        [(k-1)*n*ones(2,e);(k-2)*n*ones(2,e)] ...  
        [(k-1)*n*ones(3,e);(k-2)*n*ones(1,e)]] ) ];
    % BUG TO FIX HERE
    BDY= [BDY, ...
          bdy1 + [(k-1)*n*ones(2,bl); (k-2)*n*ones(1,bl)], ...
          bdy2 + [(k-2)*n*ones(2,bl); (k-1)*n*ones(1,bl)] ];
  end %for k

  % Now add top and bottom boundary
  BDY= [elem0, BDY, elem0+n*(ln-1) ];

  % elec_nodes is all nodes for all layers
  elec_nodes= ones(ln,1) * elec_nodes0 + ...
              (0:ln-1)'  * n*ones(1, length(elec_nodes0) );


%param.electrode = mk_electrodes( elec_nodes );
% Create the electrode structure from elec_nodes
% Currently implements point electrodes with 
%   contact impedance of near zero.
function elec_struct = mk_electrodes( elec_nodes)
   for i= 1:length( elec_nodes )
      elec_struct(i).nodes     = elec_nodes(i);
      elec_struct(i).z_contact = 0.001; % corresponds to 1 ohm
   end
   % Need to do it this way to be compatible accross versions
   if ~exist('elec_struct');
       elec_struct= [];
   end

