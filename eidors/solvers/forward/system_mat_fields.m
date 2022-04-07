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
   if length(unique(fwd_model.elems(:)))~=n
       warning('EIDORS:unused_nodes','Number of nodes inconsistent with Elements. Consider using "remove_unused_nodes"');
   end

   FFjidx= floor([0:d0*e-1]'/d0)*d1*ones(1,d1) + ones(d0*e,1)*(1:d1);
   FFiidx= [1:d0*e]'*ones(1,d1);
   FFdata = assemble_elements(d1,d0,p);

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
             compl_elec_mdl(fmdl,pp)
   d0= pp.n_dims;
   FFdata= zeros(0,d0);
   FFiidx= zeros(0,d0);
   FFjidx= zeros(0,d0);
   CCdata= zeros(0,d0);
   CCiidx= zeros(0,d0);
   CCjidx= zeros(0,d0);
  
   sidx= d0*pp.n_elem;
   cidx= (d0+1)*pp.n_elem;
   i_cem = 0; % Index into electrodes that are CEM (but not pt)
   pp.cem_boundary = pp.boundary;
   if isfield(fmdl,'system_mat_fields')
      pp.cem_boundary = [pp.cem_boundary;
          fmdl.system_mat_fields.CEM_boundary];
      pp.cem_boundary = unique(pp.cem_boundary, ...
          'rows','stable');
          % TODO, how does this get duplicated
   end
   for i= 1:pp.n_elec
      eleci = fmdl.electrode(i);
      % contact impedance zc is in [Ohm.m] for 2D or [Ohm.m^2] for 3D
      zc=  eleci.z_contact;
      [i_cem, sidx, cidx, blk]= compl_elec_mdl_i( ...
         eleci, zc, i_cem, sidx, cidx, fmdl, pp);
      FFdata= [FFdata;blk.FFdata];
      FFiidx= [FFiidx;blk.FFiidx];
      FFjidx= [FFjidx;blk.FFjidx];
      CCdata= [CCdata;blk.CCdata];
      CCiidx= [CCiidx;blk.CCiidx];
      CCjidx= [CCjidx;blk.CCjidx];
   end

function [i_cem, sidx, cidx, blk]=compl_elec_mdl_i( ...
       eleci, zc, i_cem, sidx, cidx, fmdl, pp);
   d0= pp.n_dims;
   FFd_block= sqrtm( ( ones(d0) + eye(d0) )/6/(d0-1) ); % 6 in 2D, 12 in 3D 
   FFi_block= ones(d0,1)*(1:d0);
   blk.FFdata= zeros(0,d0);
   blk.FFiidx= zeros(0,d0);
   blk.FFjidx= zeros(0,d0);
   blk.CCdata= zeros(0,d0);
   blk.CCiidx= zeros(0,d0);
   blk.CCjidx= zeros(0,d0);


   if isfield(eleci,'faces')  && ~isempty(eleci.faces)
      if ~isempty(eleci.nodes)
         eidors_msg('Warning: electrode %d has both faces and nodes',i);
      end
      bdy_nds= eleci.faces;
      [~, bdy_area] = find_electrode_bdy( ...
          bdy_nds, fmdl.nodes,[]);
   else
      [bdy_idx, bdy_area] = find_electrode_bdy( ...
        pp.cem_boundary, fmdl.nodes, eleci.nodes);
          % bdy_area is in [m] for 2D or [m^2] for 3D

      % For pt elec model, bdy_idx = [], so this doesn't run
      if isempty(bdy_idx); return; end % no action for pt elecs

      bdy_nds= pp.cem_boundary(bdy_idx,:);
   end
   i_cem = i_cem + 1;
      % 3D: [m^2]/[Ohm.m^2] = [S]
      % 2D: [m]  /[Ohm.m]   = [S]

%    for j= 1:size(bdy_nds,1);
%       blk.FFdata= [blk.FFdata; FFd_block * sqrt(bdy_area(j)/zc)];
%       blk.FFiidx= [blk.FFiidx; FFi_block' + sidx];
%       blk.FFjidx= [blk.FFjidx; FFi_block  + cidx];
% 
%       blk.CCiidx= [blk.CCiidx; FFi_block(1:2,:) + cidx];
%       blk.CCjidx= [blk.CCjidx; bdy_nds(j,:) ; (pp.n_node+i_cem)*ones(1,d0)];
%       blk.CCdata= [blk.CCdata; [1;-1]*ones(1,d0)];
%       sidx = sidx + d0;
%       cidx = cidx + d0;
%    end

   N = size(bdy_nds,1);
   blk.FFdata = kron(sqrt(bdy_area(:)/zc), FFd_block);
   blk.FFiidx = repmat(FFi_block',N,1) + kron(d0*(0:(N-1))',ones(d0)) + sidx;
   blk.FFjidx = repmat(FFi_block ,N,1) + kron(d0*(0:(N-1))',ones(d0)) + cidx;
   blk.CCiidx = repmat(FFi_block(1:2,:),N,1) + kron(d0*(0:(N-1))',ones(2,d0)) + cidx;
   blk.CCjidx = (pp.n_node+i_cem) * ones(2*N,d0);
   blk.CCjidx(1:2:end,:) = bdy_nds;
   blk.CCdata = ones(2*N, d0);
   blk.CCdata(2:2:end,:) = -1;
   sidx = sidx + N*d0;
   cidx = cidx + N*d0;
   

function  FFdata = assemble_elements(d1,d0,p);
   dfact = (d0-1)*d0;
   if 0 % OLD CODE
       FFdata= zeros(d0*p.n_elem,d1);
       for j=1:p.n_elem;
         a=  inv([ ones(d1,1), p.NODE( :, p.ELEM(:,j) )' ]);
         idx= d0*(j-1)+1 : d0*j;
         FFdata(idx,1:d1)= a(2:d1,:)/ sqrt(dfact*abs(det(a)));
       end %for j=1:ELEMs 
   end
   M3 = ones(d1,d1,p.n_elem);
   M3(2:end,:,:)= reshape(p.NODE(:,p.ELEM),d0,d1,[]);
% To stabilize inverse, remove average
% This is a multiple of row1 = 1
% This affects row one only
   M3(2:end,:,:) = bsxfun(@minus, M3(2:end,:,:), ...
              mean(M3(2:end,:,:),2));
   switch d1
      case 3; [I,D] = vectorize_3x3inv(M3);
      case 4; [I,D] = vectorize_4x4inv(M3);
      otherwise;
         error('problem dimension not understood')
   end
   Is = bsxfun(@times, I, sqrt(abs(D)/dfact));
   Is = permute(Is(:,2:d1,:),[2,3,1]);
   FFdata= reshape(Is,[],d1);

% Not vectorized yet, this is a placeholder
function [I,D] = vectorize_3x3inv(M)
   a= M(1,1,:); b= M(1,2,:); c= M(1,3,:);
   d= M(2,1,:); e= M(2,2,:); f= M(2,3,:);
   g= M(3,1,:); h= M(3,2,:); i= M(3,3,:);

   D= a.*(e.*i-f.*h) - b.*(d.*i-f.*g) + c.*(d.*h-e.*g);
   if any(abs(D) < eps) 
      warning('Determinant close to zero');
   end

   I= (1./D).* ...
      [e.*i - f.*h, c.*h - b.*i, b.*f - c.*e;
       f.*g - d.*i, a.*i - c.*g, c.*d - a.*f;
       d.*h - e.*g, b.*g - a.*h, a.*e - b.*d];
   

function [I,D] = vectorize_4x4inv(M)
% Adapted from the Mesa3D implementation of
% gluInvertMatrix(const double m[16], double invOut[16])
% Available in Mesa3D (License in MIT)
%
% positions of 1 removed in rev6202


% Precalculate pieces to speed up
   M21 = M(2,1,:); M31 = M(3,1,:); M41 = M(4,1,:);
   M22 = M(2,2,:); M32 = M(3,2,:); M42 = M(4,2,:);
   M23 = M(2,3,:); M33 = M(3,3,:); M43 = M(4,3,:);
   M24 = M(2,4,:); M34 = M(3,4,:); M44 = M(4,4,:);

   I11= M22.*M33.*M44 - ...
             M22.*M43.*M34 - ...
             M23.*M32.*M44 + ...
             M23.*M42.*M34 + ...
             M24.*M32.*M43 - ...
             M24.*M42.*M33;

   I12=-M33.*M44 + ...
             M43.*M34 + ...
             M32.*M44 - ...
             M42.*M34 - ...
             M32.*M43 + ...
             M42.*M33;

   I13= M23.*M44 - ...
             M43.*M24 - ...
             M22.*M44 + ...
             M42.*M24 + ...
             M22.*M43 - ...
             M42.*M23;

   I14=-M23.*M34 + ...
             M33.*M24 + ...
             M22.*M34 - ...
             M32.*M24 - ...
             M22.*M33 + ...
             M32.*M23;

   I21=-M21.*M33.*M44 + ...
             M21.*M43.*M34 + ...
             M23.*M31.*M44 - ...
             M23.*M41.*M34 - ...
             M24.*M31.*M43 + ...
             M24.*M41.*M33;

   I22= M33.*M44 - ...
             M43.*M34 - ...
             M31.*M44 + ...
             M41.*M34 + ...
             M31.*M43 - ...
             M41.*M33;

   I23=-M23.*M44 + ...
             M43.*M24 + ...
             M21.*M44 - ...
             M41.*M24 - ...
             M21.*M43 + ...
             M41.*M23;

   I24= M23.*M34 - ...
             M33.*M24 - ...
             M21.*M34 + ...
             M31.*M24 + ...
             M21.*M33 - ...
             M31.*M23;

   I31= M21.*M32.*M44 - ...
             M21.*M42.*M34 - ...
             M22.*M31.*M44 + ...
             M22.*M41.*M34 + ...
             M24.*M31.*M42 - ...
             M24.*M41.*M32;

   I32=-M32.*M44 + ...
             M42.*M34 + ...
             M31.*M44 - ...
             M41.*M34 - ...
             M31.*M42 + ...
             M41.*M32;

   I33= M22.*M44 - ...
             M42.*M24 - ...
             M21.*M44 + ...
             M41.*M24 + ...
             M21.*M42 - ...
             M41.*M22;

   I34=-M22.*M34 + ...
             M32.*M24 + ...
             M21.*M34 - ...
             M31.*M24 - ...
             M21.*M32 + ...
             M31.*M22;

   I41=-M21.*M32.*M43 + ...
             M21.*M42.*M33 + ...
             M22.*M31.*M43 - ...
             M22.*M41.*M33 - ...
             M23.*M31.*M42 + ...
             M23.*M41.*M32;

   I42= M32.*M43 - ...
             M42.*M33 - ...
             M31.*M43 + ...
             M41.*M33 + ...
             M31.*M42 - ...
             M41.*M32;

   I43=-M22.*M43 + ...
             M42.*M23 + ...
             M21.*M43 - ...
             M41.*M23 - ...
             M21.*M42 + ...
             M41.*M22;

   I44= M22.*M33 - ...
             M32.*M23 - ...
             M21.*M33 + ...
             M31.*M23 + ...
             M21.*M32 - ...
             M31.*M22;

   D =      I11 + ...
       M21.*I12 + ...
       M31.*I13 + ...
       M41.*I14;

   if any(abs(D) < eps) 
      warning('Determinant close to zero');
   end

   I = bsxfun(@times,...
       [I11, I12, I13, I14;
        I22, I22, I23, I24;
        I32, I32, I33, I34;
        I42, I42, I43, I44], (1./D));


function do_unit_test
   
   do_unit_test_vectorized_inv
   do_unit_test_2d_test
   do_unit_test_3d_test

function do_unit_test_vectorized_inv
   M=reshape((1:16).^(0.5),4,4); M(1,:) = 1;
   M(2:end,:,:) = M(2:end,:,:) - ...
             mean(M(2:end,:,:),2);
   [I,D] = vectorize_4x4inv(M);
   iM = inv(M);
   unit_test_cmp('vectorize_4x4inv', ...
       I(:,2:end) ,iM(:,2:end),1e-09);
   unit_test_cmp('vectorize_4x4inv',D,det(M),1e-12);
   M=reshape((1:9).^(0.5),3,3); M(1,:) = 1;
   [I,D] = vectorize_3x3inv(M);
   unit_test_cmp('vectorize_3x3inv',I,inv(M),1e-12);
   unit_test_cmp('vectorize_3x3inv',D,det(M),1e-12);

function do_unit_test_3d_test
   imdl= mk_common_model('n3r2',[16,2]); 
   fmdl= imdl.fwd_model;
   FC = system_mat_fields( fmdl );
   SZ(2) = num_nodes(fmdl) + num_elecs(fmdl);
   SZ(1) = 3*(num_elems(fmdl) + num_elecs(fmdl)*2);
     % Two faces per electrode
   unit_test_cmp('sys_mat3d', size(FC), SZ);

   E14 = [...
   0.329216512695190  -0.329216512695190  0                   0;
   0.032424996343701  -0.032424996343701 -0.506252451622855   0.506252451622855;
  -0.098764953808557                   0  0.098764953808557                   0;
   0.329216512695190  -0.329216512695190  0                   0];
   unit_test_cmp('sys_mat3dx1', FC(1:4,[1,33,64,65]), E14,1e-14);
   

function do_unit_test_2d_test

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
   Nelec = 4;
   fmdl= getfield( mk_common_model('a2C2',Nelec), 'fwd_model');

   FC = system_mat_fields( fmdl );
   unit_test_cmp('sys_mat-b2cCEM-0', Nelec, num_elecs(fmdl));
   fmdl.electrode(3).nodes = []; % no longer a node elec
   fmdl.electrode(3).faces = [1,2;2,3;3,1];
   FC = system_mat_fields( fmdl);
%spy(FC)
%for i=129:148; disp([i,find(FC(i,:))]); end
%full(FC(129:end,1:3))  


   unit_test_cmp('sys_mat-b2cCEM-1', size(FC), ...
       [num_elems(fmdl)*(size(fmdl.elems,2)-1) + 3*2 + 2*3, ...
        num_nodes(fmdl)+num_elecs(fmdl)]);

   fmdl= getfield( mk_common_model('a2C2',Nelec), 'fwd_model');
   fmdl.electrode(5) = struct('nodes',1:3,'z_contact',.01);
   fmdl.system_mat_fields.CEM_boundary = [1,2;2,3;3,1];
   FC = system_mat_fields( fmdl );

   d1=num_elems(fmdl)*(size(fmdl.elems,2)-1) + ...
        4*2 + 2*3;
   d2= num_nodes(fmdl)+num_elecs(fmdl);
   unit_test_cmp('sys_mat-b2cCEM-1', size(FC),[d1,d2]);
   unit_test_cmp('sys_mat-b2cCEM-3', FC(129:136,42:45), ...
             -13.967473716321374*kron(eye(4),[1;1]),1e-12)
   unit_test_cmp('sys_mat-ctrCEM-4', FC([137:138,141:142],end), ...
             -3.535533905932737, 1e-12);
   unit_test_cmp('sys_mat-ctrCEM-5', FC([139:140],end), ...
             -4.204482076268572, 1e-12);
