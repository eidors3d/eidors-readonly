function inv_mdl= mk_common_gridmdl( str, RM)
% MK_COMMON_MODEL: make EIT on reconstruction grids (GREIT)
%
% Usage
%      inv_mdl = mk_common_gridmdl( mdl_string, RecMtx )
%
% 2D models
%   mk_common_gridmdl('b2c', RM)  - 32x32, 16 elec
%   mk_common_gridmdl('b2d', RM)  - 32x32 - GBP shape, 16 elec

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_common_gridmdl.m,v 1.1 2008-06-11 14:52:29 aadler Exp $

if strcmp(str,'backproj')
   str= 'b2d';
   load ReconstrMatrix_GBP;
   RM= -ReconstrMatrix';
end

switch str
   case 'b2c'
      space = linspace( -1, 1, 32+1 );
      fmdl= mk_grid_model( [], space, space);
      space_avg = conv2(space, [1,1]/2,'valid'); % average of each box
      [x,y] = meshgrid( space_avg, space_avg);
      inside=  (x(:).^2 + y(:).^2)<1.15 ;
      ff_inside = find(inside); 
      map = zeros(length(space_avg)^2,1);
      map(ff_inside) = ff_inside;

   case 'b2d'
      space = linspace( -1, 1, 32+1 );
      fmdl= mk_grid_model( [], space, space);
      space_avg = conv2(space, [1,1]/2,'valid'); % average of each box
      [x,y] = meshgrid( space_avg, space_avg);
      inside=  (abs(x(:)) + abs(y(:)))<1.55 ;
      ff = find(~inside);
      fmdl.elems([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
      fmdl.coarse2fine(:,ff)= [];

      RM= RM(:,twistit); % Nodes are twisted according the goettingen format
      fmdl.nodes = -fmdl.nodes;

   otherwise
      error(['mdl_string ',str,' not understood']);
end

inv_mdl = eidors_obj('inv_model',['mk_common_gridmdl: ',str]);
inv_mdl.reconst_type= 'difference';
inv_mdl.fwd_model= fmdl;
inv_mdl.solve_use_matrix.RM = RM;
if exist('map','var')
   inv_mdl.solve_use_matrix.map = map;
end
inv_mdl.solve = @solve_use_matrix;
inv_mdl.fwd_model.normalize_measurements= 1;
[st, els]= mk_stim_patterns(16, 1, '{ad}','{ad}', {}, 10);
inv_mdl.fwd_model.meas_select= els;

function twist= twistit;
   twist= [               0+(1:13), ...
                         13+(1:13), ...
           39-(0:-1:0),  26+(1:12), ...
           52-(1:-1:0),  39+(1:11), ...
           65-(2:-1:0),  52+(1:10), ...
           78-(3:-1:0),  65+(1: 9), ...
           91-(4:-1:0),  78+(1: 8), ...
          104-(5:-1:0),  91+(1: 7), ...
          117-(6:-1:0), 104+(1: 6), ...
          130-(7:-1:0), 117+(1: 5), ...
          143-(8:-1:0), 130+(1: 4), ...
          156-(9:-1:0), 143+(1: 3), ...
          169-(10:-1:0),156+(1: 2), ...
          182-(11:-1:0),169+(1: 1), ...
          195-(12:-1:0), ...
          208-(12:-1:0) ];

