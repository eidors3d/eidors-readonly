function RM= calc_TSVD_RM(mdl,hp)
% CALC_TSVD_RM: Calculated truncated Jacobian SVD reconstruction matrix
%   RM= calc_TSVD_RM(mdl, hp)
% 
%   RM    => reconstruction matrix
%   mdl   => an inverse or forward model structure 
%   hp    => hyperparameter. Ratio between the last kept and the first SV
%            in percent (default: 1)
% 
% SEE ALSO: inv_solve_TSVD, solve_use_matrix


% (C) 2011 Bartlomiej Grychtol. Licenced under GPL v2 or v3
% $Id$

if isstr(mdl) && strcmp(mdl, 'UNIT_TEST'), do_unit_test, return,end



switch (mdl.type)
    case 'inv_model'
        fmdl = mdl.fwd_model;
        bkgnd_img = calc_jacobian_bkgnd(mdl);
    case 'fwd_model'
        fmdl = mdl;
        bkgnd_img = mk_image( fmdl,1) ;
    otherwise
        eidors_msg('calc_TSVD_RM: Object type not supported');
end

if nargin < 2
    hp = 1;
end

J= calc_jacobian(fmdl,bkgnd_img);
imdl.solve = @solve_use_matrix;
S = eidors_obj('get-cache',mdl,'TSVD_S');
if isempty(S)
    [U,S,V] = svd(J,'econ');
    eidors_obj('set-cache',mdl,'TSVD_S',S);
    eidors_obj('set-cache',mdl,'TSVD_V',V);
    eidors_obj('set-cache',mdl,'TSVD_U',U);
else
    U = eidors_obj('get-cache',mdl,'TSVD_U');
    V = eidors_obj('get-cache',mdl,'TSVD_V');
end
s = diag(S);
s = s./s(1);
N = find(s >= hp/100,1,'last'); %disp(N);
t= 1:N; % tsvd
RM = eidors_obj('get-cache',mdl,'TSVD_RM',N);
if isempty(RM)
    RM  = (V(:,t)/S(t,t))*U(:,t)';
    eidors_obj('set-cache',mdl,'TSVD_RM',RM,N);
end



function do_unit_test
mdl = mk_common_model('b3t2r',[16,1]);
RM = calc_TSVD_RM(mdl);