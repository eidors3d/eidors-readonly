function FT= system_mat_2p5d_fields( fwd_model )
% SYSTEM_MAT_2P5D_FIELDS: fields (elem to nodes) fraction of system mat
% FC= system_mat_fields( fwd_model )
% input: 
%   fwd_model = forward model
% output:
%   FT:        s_mat = C' * T * conduct * C = FT' * conduct * FT;
%   where FC'*S*FC + k^2*FT'*S*FT, with k=real#, forms part of the Fourier integral

% (C) 2008, 2015, 2016 Andy Adler, Alistair Boyle. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

copt.cache_obj = mk_cache_obj(fwd_model);
copt.fstr = 'system_mat_2p5d_fields';
copt.log_level = 4;
FT = eidors_cache(@calc_system_mat_2p5d_fields,{fwd_model},copt );

% only cache stuff which is really relevant here
function cache_obj = mk_cache_obj(fwd_model)
   cache_obj.elems       = fwd_model.elems;
   cache_obj.nodes       = fwd_model.nodes;
   try
   cache_obj.electrode   = fwd_model.electrode; % if we have it
   end
   cache_obj.type        = 'fwd_model';
   cache_obj.name        = ''; % it has to have one

function FT = calc_system_mat_2p5d_fields( fwd_model )
   p= fwd_model_parameters( fwd_model, 'skip_VOLUME' );
   d0= p.n_dims+0;
   d1= p.n_dims+1;
   e= p.n_elem;
   n= p.n_node;
   assert(d0 == 2, 'only supports 2D meshes for 2.5D system matrices');

   dfact = (d0-1)*d0;
   SSe2 = sqrtm((ones(3) + eye(3)) / 12);
   area=zeros(e,1);
   for j=1:e
      a = inv([ ones(d1,1), p.NODE( :, p.ELEM(:,j) )' ]);
      area(j) = 1 / sqrt(dfact*abs(det(a)));
   end %for j=1:ELEMs

   % NOTE: CEM electrodes are taken care of by the system_mat_fields produced
   % FC matrix where we sum FC'*S*FC + k^2*FT'*S*FT for 2.5D solutions, so we
   % just need to have the correctly sized FT matrix, so that the FC and FT
   % products may be added directly
   n_cem = 0; % count CEM electrodes, add 1 column & row per CEM (all zeros)
   for i=1:length(fwd_model.electrode)
      n_cem = n_cem + (length(fwd_model.electrode(i).nodes) > 1);
   end
   [FFiidx,FFjidx,FFdata] = find(kron(spdiag(area),SSe2));
   FF= sparse(FFiidx(:), FFjidx(:), FFdata(:));
   CC= sparse((1:d1*e)', p.ELEM(:), ones(d1*e,1), d1*e, n+n_cem);

   FT= FF*CC;

function do_unit_test
   imdl=  mk_common_model('a2c2',16);
   FT = system_mat_2p5d_fields( imdl.fwd_model);
   unit_test_cmp('sys_mat1 pem', size(FT), [192,41]);
   ss = 58.787753826797278;
   unit_test_cmp('sys_mat2 pem', FT(1:2,:), [[4,1,1,0;1,4,1,0], zeros(2,37)]/ss, 1e-14);

   imdl=  mk_geophysics_model('h2a', 16);
   FC = system_mat_fields( imdl.fwd_model);
   FT = system_mat_2p5d_fields( imdl.fwd_model);
   unit_test_cmp('sys_mat1 cem', size(FT, 2), size(FC, 2));
   unit_test_cmp('sys_mat2 cem', FT(:, end-15:end), FT(:, end-15:end)*0);
