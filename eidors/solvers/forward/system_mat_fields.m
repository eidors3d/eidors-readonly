function FC= system_mat_fields( fwd_model )
% SYSTEM_MAT_FIELDS: fields (elem to nodes) fraction of system mat
% FC= system_mat_fields( fwd_model )
% input: 
%   fwd_model = forward model
% output:
%   FC:        s_mat= C' * S * conduct * C = FC' * conduct * FC;
%
% system_mat_fields detects whether an electrode is a CEM by checking
%  1) whether it has more than 1 node, and 2) if it's on the boundary
% If you want an internal CEM that's not on the boundary, then it
% has to be specified by using
%    fmdl.system_mat_fields.CEM_boundary = [extra boundary ]

% (C) 2008-2018 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

copt.cache_obj = mk_cache_obj(fwd_model);
copt.fstr = 'system_mat_fields';
copt.log_level = 4;
FC= eidors_cache(@calc_system_mat_fields,{fwd_model},copt );


% only cache stuff which is really relevant here
function cache_obj = mk_cache_obj(fwd_model)
   cache_obj.elems       = fwd_model.elems;
   cache_obj.nodes       = fwd_model.nodes;
   try
   cache_obj.electrode   = fwd_model.electrode; % if we have it
   end
   cache_obj.type        = 'fwd_model';
   cache_obj.name        = ''; % it has to have one

function FC= calc_system_mat_fields( fwd_model )
   p= fwd_model_parameters( fwd_model, 'skip_VOLUME' );
   d0= p.n_dims+0;
   d1= p.n_dims+1;
   e= p.n_elem;
   n= p.n_node;

   FFjidx= floor([0:d0*e-1]'/d0)*d1*ones(1,d1) + ones(d0*e,1)*(1:d1);
   FFiidx= [1:d0*e]'*ones(1,d1);
   FFdata= zeros(d0*e,d1);
   dfact = (d0-1)*d0;
   for j=1:e
     a=  inv([ ones(d1,1), p.NODE( :, p.ELEM(:,j) )' ]);
     idx= d0*(j-1)+1 : d0*j;
     FFdata(idx,1:d1)= a(2:d1,:)/ sqrt(dfact*abs(det(a)));
   end %for j=1:ELEMs 

if 0 % Not complete electrode model
   FF= sparse(FFiidx,FFjidx,FFdata);
   CC= sparse((1:d1*e),p.ELEM(:),ones(d1*e,1), d1*e, n);
else
   [F2data,F2iidx,F2jidx, C2data,C2iidx,C2jidx] = ...
             compl_elec_mdl(fwd_model,p);
   FF= sparse([FFiidx(:); F2iidx(:)],...
              [FFjidx(:); F2jidx(:)],...
              [FFdata(:); F2data(:)]);
   
   CC= sparse([(1:d1*e)';    C2iidx(:)], ...
              [p.ELEM(:);   C2jidx(:)], ...
              [ones(d1*e,1); C2data(:)]);
end

FC= FF*CC;

% Add parts for complete electrode model
function [FFdata,FFiidx,FFjidx, CCdata,CCiidx,CCjidx] = ...
             compl_elec_mdl(fwd_model,pp)
   d0= pp.n_dims;
   FFdata= zeros(0,d0);
   FFd_block= sqrtm( ( ones(d0) + eye(d0) )/6/(d0-1) ); % 6 in 2D, 12 in 3D 
   FFiidx= zeros(0,d0);
   FFjidx= zeros(0,d0);
   FFi_block= ones(d0,1)*(1:d0);
   CCdata= zeros(0,d0);
   CCiidx= zeros(0,d0);
   CCjidx= zeros(0,d0);
  
   sidx= d0*pp.n_elem;
   cidx= (d0+1)*pp.n_elem;
   i_cem = 0; % Index into electrodes that are CEM (but not pt)
   cem_boundary = pp.boundary;
   if isfield(fwd_model,'system_mat_fields')
      cem_boundary = [cem_boundary;fwd_model.system_mat_fields.CEM_boundary];
   end
   for i= 1:pp.n_elec
      eleci = fwd_model.electrode(i);
      % contact impedance zc is in [Ohm.m] for 2D or [Ohm.m^2] for 3D
      zc=  eleci.z_contact;
      [bdy_idx, bdy_area] = find_electrode_bdy( ...
          cem_boundary, fwd_model.nodes, eleci.nodes );
          % bdy_area is in [m] for 2D or [m^2] for 3D

      % For pt elec model, bdy_idx = [], so this doesn't run
      if isempty(bdy_idx); continue; end % no action for pt elecs
      i_cem = i_cem + 1;

      for j= 1:length(bdy_idx);
         bdy_nds= cem_boundary(bdy_idx(j),:);

         % 3D: [m^2]/[Ohm.m^2] = [S]
         % 2D: [m]  /[Ohm.m]   = [S]
         FFdata= [FFdata; FFd_block * sqrt(bdy_area(j)/zc)];
         FFiidx= [FFiidx; FFi_block' + sidx];
         FFjidx= [FFjidx; FFi_block  + cidx];

         CCiidx= [CCiidx; FFi_block(1:2,:) + cidx];
         CCjidx= [CCjidx; bdy_nds ; (pp.n_node+i_cem)*ones(1,d0)];
         CCdata= [CCdata; [1;-1]*ones(1,d0)];
         sidx = sidx + d0;
         cidx = cidx + d0;
      end
      
   end

function do_unit_test
   imdl=  mk_common_model('a2c2',16);
   FC = system_mat_fields( imdl.fwd_model);
   unit_test_cmp('sys_mat1', size(FC), [128,41]);
   unit_test_cmp('sys_mat2', FC(1:2,:), [[0,-1,1,0;-2,1,1,0], zeros(2,37)]/2, 1e-14);

   % THis is a 45 degree rotation of the previous
   imdl=  mk_common_model('a2c0',16);
   FC2= system_mat_fields( imdl.fwd_model);
   M = sqrt(.5)*[1,-1;1,1];
   unit_test_cmp('sys_mat3', M*FC2(1:2,:), [[0,-1,1,0;-2,1,1,0], zeros(2,37)]/2, 1e-14);

   % Check we can handle partial CEM/pt electrodes
   imdl = mk_common_model('b2s',4); fmdl = imdl.fwd_model;
   fmdl.electrode([2,4])= []; % remove to simplify
   fmdl.stimulation = stim_meas_list([1,2,1,2]);

   FF = system_mat_fields(fmdl);
   unit_test_cmp('sys_mat_b2s_1a', size(FF), [256,81]);
   vh = fwd_solve( mk_image(fmdl,1) ); 
   unit_test_cmp('sys_mat_b2s_1b', vh.meas, 2.182865407049302, 1e-12);

   fmdl.electrode(2).nodes = [4:6];
   FF = system_mat_fields(fmdl);
   unit_test_cmp('sys_mat_b2s_2a', size(FF), [260,82]);
   vh = fwd_solve( mk_image(fmdl,1) ); 
   unit_test_cmp('sys_mat_b2s_2b', vh.meas, 1.807627806615849, 1e-12);

   fmdl.electrode(1).nodes = [76:78];
   FF = system_mat_fields(fmdl);
   unit_test_cmp('sys_mat_b2s_3a', size(FF), [264,83]);
   vh = fwd_solve( mk_image(fmdl,1) ); 
   unit_test_cmp('sys_mat_b2s_3b', vh.meas, 1.432226638247073, 1e-12);

   imdl=  mk_common_model('a2C2',4); fmdl=imdl.fwd_model;
   fmdl.nodes(1,:) = [];
   fmdl.gnd_node = 2;
   fmdl.elems(1:4,:) = [];
   fmdl.elems = fmdl.elems - 1;
   fmdl.boundary = find_boundary(fmdl);
   FC = system_mat_fields( fmdl);
   unit_test_cmp('sys_mat-bdyCEM-1', size(FC),[128,44]);
   unit_test_cmp('sys_mat-bdyCEM-2', FC(121:end,41:end), ...
             -13.967473716321374*kron(eye(4),[1;1]),1e-12)

% Add a CEM via a boundary
   fmdl.electrode(5) = struct('nodes',1:4,'z_contact',.01);
   FC = system_mat_fields( fmdl);
   unit_test_cmp('sys_mat-ctrCEM-1', size(FC),[136,45]);
   unit_test_cmp('sys_mat-ctrCEM-2', FC(121:128,41:44), ...
             -13.967473716321374*kron(eye(4),[1;1]),1e-12)
   unit_test_cmp('sys_mat-ctrCEM-2', FC(129:end,end), ...
             -4.204482076268572,1e-12)

% Add a CEM via a boundary
   imdl=  mk_common_model('a2C2',4); fmdl = imdl.fwd_model;
   fmdl.electrode(5) = struct('nodes',1:3,'z_contact',.01);
   fmdl.system_mat_fields.CEM_boundary = [1,2;2,3;3,1];
   FC = system_mat_fields( fmdl );
   unit_test_cmp('sys_mat-b2cCEM-1', size(FC),[142,46]);
   unit_test_cmp('sys_mat-b2cCEM-2', FC(129:136,42:45), ...
             -13.967473716321374*kron(eye(4),[1;1]),1e-12)
   unit_test_cmp('sys_mat-ctrCEM-3', FC([137:138,141:142],end), ...
             -3.535533905932737, 1e-12);
   unit_test_cmp('sys_mat-ctrCEM-4', FC([139:140],end), ...
             -4.204482076268572, 1e-12);
