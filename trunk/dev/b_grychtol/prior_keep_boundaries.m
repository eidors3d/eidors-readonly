function Reg = prior_keep_boundaries(inv_model)

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

switch inv_model.type
  case 'inv_model'; fwd_model = inv_model.fwd_model;
  case 'fwd_model'; fwd_model = inv_model;
  otherwise; error('PRIOR_KEEP_BOUNDARIES requires input type of inv_model or fwd_model');
end

try 
    func = inv_model.prior_keep_boundaries.func;
catch
    func = eidors_default('get','calc_RtR_prior');
    eidors_msg('@@@ using eidors_default prior function');
end

fmdl = fwd_model;
R = []; C =[]; V = [];
for i = 1:length(fwd_model.mat_idx)
    m = fwd_model.mat_idx{i};
    fmdl.elems = fwd_model.elems(m,:);
    P = feval(func, fmdl);
    [r c v] = find(P);
    R = [R; m(r(:))];
    C = [C; m(c(:))];
    V = [V; v(:)];
end
nel = length(fwd_model.elems);
Reg = sparse(R,C,V,nel,nel);

function do_unit_test
mdl = mk_library_model('pig_23kg_16el_lungs');
[stim msel] = mk_stim_patterns(16,1,4,'{ad}',{},1);
mdl.stimulation = stim;
mdl.meas_select = msel;
img = mk_image(mdl,1);
img.elem_data(mdl.mat_idx{2}) = 0.3;
vh = fwd_solve(img);
img.elem_data(mdl.mat_idx{2}) = 0.2;
vi = fwd_solve(img);
imdl = select_imdl(mdl, {'Basic GN dif'});
% imdl = add_rec_model(imdl, [128 128]);
% imdl.prior_use_fwd_not_rec = 1;
% imdl.fwd_model = rmfield(imdl.fwd_model,'coarse2fine');
rimg = inv_solve(imdl, vi, vh);
subplot(121)
show_slices(rimg,[inf inf 0.5]);
imdl.RtR_prior = @prior_keep_boundaries;
subplot(122)
rimg = inv_solve(imdl, vi, vh);
show_slices(rimg,[inf inf 0.5]);

% R = prior_keep_boundaries(mdl);