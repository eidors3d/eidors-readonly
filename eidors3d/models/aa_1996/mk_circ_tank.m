function param= mk_circ_tank(rings, levels, n_elec, n_planes )
%MK_CIRC_TANK: make a cylindrical tank FEM geometry in 2D or 3D
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
% 
% rings:  number of horizontal plane rings (divisible by 4)
% levels: vector of vertical placement of levels
%     for 2D mesh, levels = []
% n_elec: number of electrodes in each horiz plane (divisible by 4)
% n_planes: number of planes of electrodes (divisible by 4)
%
% mk_circ_tank creates simple, point electrodes. If you wish
%   to implement more complete models, functions will need to
%   be overridden in this model
%
% output:
%  param.name        Model name (if known) 
%  param.nodes       position of FEM nodes (Nodes x Dims) 
%  param.elems       definition of FEM elements (Elems x Dims+1) 
%  param.boundary    nodes of element faces on the medium surface 
%  param.gnd_node    Number of node connected to ground 
%  param.electrode   Vector (Num_elecs x 1) of electrode models (elec_model) 

if rem(rings,4) ~= 0 || rem(n_elec,4) ~= 0;
   error('parameter rings and n_elec must be divisible by 4');
end


[elem, node, bdy, elec_nodes] = mk_2D_model( rings, n_elec);
if ~isempty( levels )
   [elem, node, bdy, elec_nodes] = mk_3D_model( elem, node, levels, bdy );
end

param.name= sprintf('EIT FEM by mk_circ_tank with N=%d levs=%d', ...
                    rings, length(levels) );
param.nodes = node';
param.elems = elem';
param.boundary = bdy';
param.gnd_node = 1; % node at bottom and center of the tank
%param.electrode = mk_electrodes( elec_nodes );


% Create a simple 2D regular mesh, based on N circular rings
%   and n_elec electrodes
function [ELEM, NODE, bdy_nodes, elec_nodes] = mk_2D_model( N, n_elec );
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
  elec_nodes= ext+N*4/n_elec*([0:n_elec-1]);
 

% 'extrude' a 2D model defined by ELEM and NODE into a 3D model
% levels are defined by 'niveaux', 
% 2D parameters are ELEM, NODE, and bdy
function [ELEM, NODE, BDY, elec_nodes] = mk_3D_model(elem0, node0, niveaux, bdy );
  d= size(elem0,1);       %dimentions+1
  n= size(node0,2);       %NODEs
  e= size(elem0,2);       %ELEMents     

  elem= [elem0([1 1 2 3],:), ...
         elem0([2 1 2 3],:), ...
         elem0([3 2 1 3],:)]; 
  NODE= [node0; niveaux(1)*ones(1,n) ];
  ELEM= [];
  bdy1= [bdy;bdy(1,:)];
  bdy2= [bdy;bdy(2,:)];
  bl= size(bdy,2);
  BDY = [];
 
  ln= length(niveaux);
  for k=2:ln
    NODE=[NODE  [node0; niveaux(k)*ones(1,n)] ];
    BDY= [BDY, ...
          bdy1 + [(k-1)*n*ones(2,bl); (k-2)*n*ones(1,bl)], ...
          bdy2 + [(k-2)*n*ones(2,bl); (k-1)*n*ones(1,bl)] ];
    ELEM= [ELEM (elem + ...
       [[(k-1)*n*ones(1,e);(k-2)*n*ones(3,e)] ...
        [(k-1)*n*ones(2,e);(k-2)*n*ones(2,e)] ...  
        [(k-1)*n*ones(3,e);(k-2)*n*ones(1,e)]] ) ];
  end %for k

  % Now add top and bottom boundary 
  BDY= [elem0, BDY, elem0+n*(ln-1) ];

  elec_nodes = [];

% MES= MES + floor(length(niveaux)/2)*n;
% cour(1:elec,:)= cour(1:elec,:)+ floor(length(niveaux)/2)*n;

%param.electrode = mk_electrodes( elec_nodes );
% Create the electrode structure from elec_nodes
% Currently implements point electrodes with 
%   contact impedance of zero.
function elec_struct = mk_electrodes( elec_nodes)
   for i= 1:length( elec_nodes )
      elec_struct(i).nodes     = elec_nodes(i);
      elec_struct(i).z_contact = 0;
   end

