function RM = unpack_reconst_matrix(packed_matrix, Nelec, Ngrid, options);
% UNPACK_RECONST_MATRIX: unpack a compressed, stored reconstruction matrix
% Reconstruction matrices, especially on a circular body, have many
%  symmetries (left-right, up-down, flip, and recoprocity = 16x).
% To save space, a matrix can be stored in this format, and then
%  used to save space. This function is used to rebuild the matrix
%
% RM = unpack_reconst_matrix(packed_matrix, Nelec, Ngrid);
%   Nelec = number of electrodes
%   Ngrid = number of grid points
%
% Example:
%  load GREIT_v10_Circ_Matrix.mat
%  RM = unpack_reconst_matrix(GREIT_v10_Circ_Matrix, 16, 32);


% (C) 2008-2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

elec_idx = 1:Nelec;
grid_idx = 1:Ngrid; grid_idx = grid_idx - mean(grid_idx);
Nelec2 = Nelec/2;
Nelec34= Nelec*3/4;

   % Take a slice
   [x,y]= meshgrid( elec_idx, elec_idx);
   ss1 = (y-x)>1 & (y-x)<15;
   sel1 = abs(x-y)>1 & abs(x-y)<15;
   
   [x,y]= meshgrid( grid_idx, grid_idx);
   ss2 = abs(x-y)<25 & abs(x+y)<25 ...
       & x<0 & y<0 & x>=y ;
   sel2 = abs(x-y)<25 & abs(x+y)<25;
 
   % Build up
   BP  = zeros(Nelec^2, Ngrid^2);
   BP(ss1,ss2) = packed_matrix;
   BP  = reshape(BP, Nelec,Nelec,Ngrid,Ngrid);

   % Reciprocity
   BP  = BP + permute(BP, [2,1,3,4]);

   % FLIP LR
   el= Nelec:-1:1;
   BP= BP + BP(el,el,[Ngrid:-1:1],:);
   % FLIP UD
   el= [Nelec2:-1:1,Nelec:-1:Nelec2+1];
   BP= BP + BP(el,el,:,[Ngrid:-1:1]);
   % Transpose
   el= [Nelec34:-1:1,Nelec:-1:Nelec34+1];
   BP= BP + permute(BP(el,el,:,:), [1,2,4,3]);

   % Final UD flip to match radiological view (upward toward patient)
   % Here electrodes are connected CW starting from TDC
   BP= BP(:,:,:,[Ngrid:-1:1]);

   RM= reshape(BP, Nelec^2, [])';
   RM= RM(:,sel1);
% This creates the diamond shape, but we want to leave shape choice later
%  RM= RM(sel2,sel1);
