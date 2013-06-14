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
fmdl = rmfield(fmdl, 'coarse2fine');

R = []; C =[]; V = [];

try 
   weights = inv_model.prior_keep_boundaries.weights;
catch
   weights = ones(length(fwd_model.mat_idx),1);
end

for i = 1:length(fwd_model.mat_idx)
    m = fwd_model.mat_idx{i};
    fmdl.elems = fwd_model.elems(m,:);
    P = weights(i)*feval(func, fmdl);
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
img1 = img;
% img.elem_data(mdl.mat_idx{2}) = 0.2;
select_fcn = inline('(x+0.2).^2+(y+0.2).^2+(z-0.5).^2<0.2^2','x','y','z');
memb_frac = elem_select( img.fwd_model, select_fcn);
img.elem_data = img.elem_data + memb_frac*0.1;
vi = fwd_solve(img);
img2 = img;
imdl = select_imdl(mdl, {'Basic GN dif'});
imdl = add_rec_model(imdl, [64 64]);
% imdl.prior_use_fwd_not_rec = 1;
% imdl.fwd_model = rmfield(imdl.fwd_model,'coarse2fine');

subplot(231)
rimg = inv_solve(imdl, vi, vh);
show_slices(rimg,[inf inf 0.5]);

subplot(232)
imdl.RtR_prior = @prior_keep_boundaries;
rimg = inv_solve(imdl, vi, vh);
show_slices(rimg,[inf inf 0.5]);

subplot(233)
imdl.prior_keep_boundaries.func = @prior_tikhonov;
imdl.prior_keep_boundaries.weights = [4 0.5];
rimg = inv_solve(imdl, vi, vh);
show_slices(rimg,[inf inf 0.5]);

subplot(234)
imdl.RtR_prior = @prior_combine;
imdl.prior_combine.prior1.func = @prior_keep_boundaries;
imdl.prior_combine.prior1.prior_keep_boundaries.func = @prior_laplace;
imdl.prior_combine.prior1.prior_keep_boundaries.weights = [0 .2];
imdl.prior_combine.prior2.func = @prior_keep_boundaries;
imdl.prior_combine.prior2.prior_keep_boundaries.func = @prior_tikhonov;
imdl.prior_combine.prior2.prior_keep_boundaries.weights = [4 .1];
rimg = inv_solve(imdl, vi, vh);
show_slices(rimg,[inf inf 0.5]);

subplot(235)
% imdl.prior_keep_boundaries.weights = [1 1];
imdl = select_imdl(imdl,{'Elec Move GN'});
rimg = inv_solve(imdl, vi, vh);
show_slices(rimg,[inf inf 0.5]);
% show_slices_move(rimg,[inf inf 0.5]);


subplot(236)
difimg = img;
difimg.elem_data = img2.elem_data - img1.elem_data;
show_slices(difimg,[inf inf 0.5]);
% imdl.prior_keep_boundaries.weights = [5 0.5];
% rimg = inv_solve(imdl, vi, vh);
% show_slices(rimg,[inf inf 0.5]);

% R = prior_keep_boundaries(mdl);