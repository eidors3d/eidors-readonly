function sol = inv_solve_TSVD(inv_model, data1, data2)
% INV_SOLVE_TSVD: inverse solver based on truncatated SVD
%   img= inv_solve_TSVD(inv_model, data1, data2)
% 
%   img        => output image (or vector of images)
%   data1      => differential data at earlier time
%   data2      => differential data at later time
%   inv_model  => inverse model struct. Requires:
%                 inv_model.hyperparameter.value OR
%                 inv_model.hyperparameter.func (and possibly others)
%
% SEE ALSO: calc_TSVD_RM, calc_hyperpameter, solve_use_matrix

% (C) 2011 Bartlomiej Grychtol. Licenced under GPL v2 or v3
% $Id$

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'), do_unit_test, return, end;

hp = calc_hyperparameter(inv_model);
inv_model.solve_use_matrix.RM = calc_TSVD_RM(inv_model, hp);
sol = solve_use_matrix(inv_model,data1,data2);


function do_unit_test
    % get a fwd_model
    mdl = mk_common_model('c3cr',16);fmdl = mdl.fwd_model; clear mdl;
    % homogenous measurement
    img  = mk_image(fmdl,1);
    vh = fwd_solve(img);
    % inhomogeneous measurement
    str = sprintf('(x-%f).^2+(y-%f).^2+(z-%f).^2<%f^2',[-0.3 0.2 0 0.1]);
    select_fcn = inline(str,'x','y','z');
    e = elem_select(img.fwd_model, select_fcn);
    ind = find(e);
    img.elem_data(ind) = img.elem_data(ind) - 0.25*e(ind);
    vi = fwd_solve(img);

    % build an inverse model
    imdl.name= 'Test TSVD model';
    imdl.type= 'inv_model';
    imdl.solve= @inv_solve_TSVD;
    imdl.hyperparameter.value = 1;
    imdl.jacobian_bkgnd.value = 1;
    imdl.reconst_type= 'difference';
    imdl.fwd_model = fmdl;
    
    % solve
    rimg = inv_solve(imdl, vh, vi);
    show_slices(rimg,[inf,inf,0]);