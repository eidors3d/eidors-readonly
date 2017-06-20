function mr = calc_model_reducton(fmdl)
% calc_model_reduction: calculate the fields for a reduced model
%    which should speed up forward calculations
% 
% mr = calc_model_reducton(fmdl)
%  where
%    mr.main_region = vector, and 
%    mr.regions = struct
%
% Model Reduction: use precomputed fields to reduce the size of
%    the forward solution. Nodes which are 1) not used in the output
%    (i.e. not electrodes) 2) all connected to the same conductivity via
%    the c2f mapping are applicable.

% see: Model Reduction for FEM Forward Solutions, Adler & Lionheart, EIT2016

% $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

switch fmdl(1).type
  case 'image';      fmdl = fmdl.fwd_model;
  case 'inv_model';  fmdl = fmdl.fwd_model;
  case 'fwd_model';  fmdl = fmdl;
  otherwise;
      error('can''t process model of type %s', fmdl.type );
end

mr.main_region = logical(ones(1:num_nodes(fmdl))'); 
mr.regions     = struct();
if ~isfield(fmdl,'coarse2fine'); return; end % not needed without


copt.cache_obj.elems = fmdl.elems;
copt.cache_obj.nodes = fmdl.nodes;
copt.cache_obj.coarse2fine = fmdl.coarse2fine;
copt.fstr = 'calc_model_reducton';
copt.log_level = 4;
mr= eidors_cache(@calc_model_reducton_fn,{fwd_model},copt );

function mr = calc_model_reducton_fn(fwd_model);

function do_unit_test
   fmdl = mk_common_model('a2c0',16); img = mk_image(fmdl,1);
   img.fwd_model.electrode = img.fwd_model.electrode([15,16,1:3]);
   stim = mk_stim_patterns(5,1,[0,1],[0,1],{},1);
   img.fwd_model.stimulation = stim;

   pos = interp_mesh(img.fwd_model);
   bl = find(pos(:,1)<=0.2 & pos(:,2)<=0.5);
   br = find(pos(:,1)>=0.2 & pos(:,2)<=0.5);
   img.elem_data(bl) = 0.9;
   img.elem_data(br) = 1.1;
   pgnode = zeros(1,mdl_dim(img)); pgnode(2) = 0.8;
   [~,img.fwd_model.gnd_node] = min(sum( ...
               bsxfun(@minus,img.fwd_model.nodes,pgnode).^2,2));
   c2f = 1+(1:num_elems(img));
   c2f(bl) = 1;
   c2f(br) = 2+num_elems(img);
   [~,idx2,idx] = unique(c2f); nc = length(idx2);
   c2f = sparse(1:num_elems(img),idx,1);

   img.elem_data = c2f*linspace(.7,1.3,length(idx2))';

   vs1 = fwd_solve( img);

   imgmr = img; imgmr.model_reduction = @calc_model_reduction;
   vs2 = fwd_solve( imgmr);
   unit_test_cmp('1-2',vs1.meas,vs2.meas,1e-10);


