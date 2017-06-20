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

mr.main_region = logical(ones(1,num_nodes(fmdl))'); 
mr.regions     = struct([]);
if ~isfield(fmdl,'coarse2fine'); return; end % not needed without


copt.cache_obj.elems = fmdl.elems;
copt.cache_obj.nodes = fmdl.nodes;
copt.cache_obj.coarse2fine = fmdl.coarse2fine;
copt.fstr = 'calc_model_reducton';
copt.log_level = 4;
mr= eidors_cache(@calc_model_reducton_fn,{fwd_model},copt );

function mr = calc_model_reducton_fn(fwd_model);
   ne = num_elems(fwd_model);
   el = fwd_model.elems; ei=(1:ne)';
   e2n = sparse(el,ei*ones(1,elem_dim(fwd_model)+1),1);

   c2n = e2n*fwd_model.coarse2fine;

% get the number of the c2n mapping. But only work on those with one
   ac2n = c2n>0;
   ac2n(sum(ac2n,2)>1,:) = 0;
   region = ac2n*(1:nc)';
   region0 = region==0;
   FC = system_mat_fields(img.fwd_model); S= FC'*FC;
   region0(end+1:size(FC,2)) = logical(1);
   k=0;
   for i=1:max(region);
      ff=find(region==i);
      if length(ff)>1;
         imgn.node_data(ff) = 1+0.5*(k)/5; k=k+1;
        %invEi = S(region0,ff)*inv(S(ff,ff))*S(ff,region0);
         invEi =(S(region0,ff)/S(ff,ff))*S(ff,region0);
         regions(k) = struct('nodes',ff,'field',i,'invE',invEi);
      end;
   end

   mr.main_region = region0;
   mr.regions = regions;

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

   imgmr = img; imgmr.fwd_model.model_reduction = @calc_model_reduction;
   vs2 = fwd_solve( imgmr);
   unit_test_cmp('1-2',vs1.meas,vs2.meas,1e-10);

   imgm3 = img; img.fwd_model.model_reduction = calc_model_reduction( fmdl);
   vs3 = fwd_solve( imgm3);
   unit_test_cmp('1-3',vs1.meas,vs3.meas,1e-10);

   imgm4 = img; img.fwd_model.model_reduction = @calc_model_reduction;
   vs4 = fwd_solve( imgm4);
   unit_test_cmp('1-4',vs1.meas,vs3.meas,1e-10);
