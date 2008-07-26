function test_c2f_jacobian
% Test calc of jacobian given coarse to fine mapping
%
% (C) Andy Adler 2008. Licenced under GPL v2 or v3
% $Id$
eidors_cache clear

test1_np_3d
test2_aa_2d
test3_aa_3d
test4_npaa_3d

function test1_np_3d
% Fine model
imdl = mk_common_model('n3r2');
f1mdl = imdl.fwd_model;

% Create 2D FEM of all NODES with z=0
n2d = f1mdl.nodes( (f1mdl.nodes(:,3) == 0), 1:2);
e2d = delaunayn(n2d);
c_mdl = eidors_obj('fwd_model','2d','elems',e2d,'nodes',n2d);

% Create f_mdl with coarse2fine mapping
c2f= mk_coarse_fine_mapping( f1mdl, c_mdl );
f2mdl= f1mdl;
f2mdl.coarse2fine = c2f;

% Just in case - define reconstruction
imdl.solve=       @np_inv_solve;
imdl.hyperparameter.value = 1e-3;
imdl.R_prior= @np_calc_image_prior;
imdl.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
imdl.jacobian_bkgnd.value= .6;
imdl.reconst_type= 'difference';
imdl.fwd_model= f1mdl; % fine model

img = calc_jacobian_bkgnd( imdl);

t=cputime;
J1= np_calc_jacobian( f1mdl, img)*c2f;
fprintf('np_calc - full: t=%f\n',  cputime-t  );
n_J1 = norm(J1,'fro');

t=cputime;
J2= np_calc_jacobian( f2mdl, img);
fprintf('np_calc - c2f: t=%f\n',  cputime-t  );

n_J1_J2 = norm(J1-J2,'fro')/ n_J1;
if n_J1_J2 > 1e-4; 
   warning('J1-J2 is too big');
end


t=cputime;
J2p= perturb_jacobian( f2mdl, img);
fprintf('perturb_j - c2f: t=%f\n',  cputime-t  );

n_J1_J2p = norm(J1-J2p,'fro')/ n_J1;
if n_J1_J2p > 1e-4; 
   warning('J1-J2p is too big');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test2_aa_2d

imdl= mk_common_model('a2c2',16);
cmdl = imdl.fwd_model;
imdl= mk_common_model('c2c0',16);
f1mdl = imdl.fwd_model;

c2f= mk_coarse_fine_mapping(f1mdl,cmdl);
img= calc_jacobian_bkgnd(imdl);

t=cputime;
J1= aa_calc_jacobian( f1mdl, img)*c2f;
fprintf('aa_calc - full: t=%f\n',  cputime-t  );
n_J1 = norm(J1,'fro');

f2mdl= f1mdl;
f2mdl.coarse2fine = c2f;

%%%%%%%%%% J2
t=cputime;
J2= aa_calc_jacobian( f2mdl, img);
fprintf('aa_calc - c2f: t=%f\n',  cputime-t  );

n_J1_J2 = norm(J1-J2,'fro')/ n_J1;
if n_J1_J2 > 1e-4; 
   warning('J1-J2 is too big');
keyboard
end


%%%%%%%%%% J2p
t=cputime;
J2p= perturb_jacobian( f2mdl, img);
fprintf('perturb_j - c2f: t=%f\n',  cputime-t  );

n_J1_J2p = norm(J1-J2p,'fro')/ n_J1;

if n_J1_J2p > 1e-4; 
   warning('J1-J2p is too big');
end

%%%%%%%%%% J1p
t=cputime;
J1p= perturb_jacobian( f1mdl, img)*c2f;
fprintf('perturb_j - full: t=%f\n',  cputime-t  );

n_J1_J1p = norm(J1-J1p,'fro')/ n_J1;
if n_J1_J1p > 1e-4; 
   warning('J1-J1p is too big');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test3_aa_3d

imdl= mk_common_model('a2c0',16);
cmdl = imdl.fwd_model;
imdl= mk_common_model('a3cr',16);
f1mdl = imdl.fwd_model;

c2f= mk_coarse_fine_mapping(f1mdl,cmdl);
img= calc_jacobian_bkgnd(imdl);

t=cputime;
J1= aa_calc_jacobian( f1mdl, img)*c2f;
fprintf('aa_calc - full: t=%f\n',  cputime-t  );
n_J1 = norm(J1,'fro');

f2mdl= f1mdl;
f2mdl.coarse2fine = c2f;

%%%%%%%%%% J2
t=cputime;
J2= aa_calc_jacobian( f2mdl, img);
fprintf('aa_calc - c2f: t=%f\n',  cputime-t  );

n_J1_J2 = norm(J1-J2,'fro')/ n_J1;
if n_J1_J2 > 1e-4; 
   warning('J1-J2 is too big');
keyboard
end


%%%%%%%%%% J2p
t=cputime;
J2p= perturb_jacobian( f2mdl, img);
fprintf('perturb_j - c2f: t=%f\n',  cputime-t  );

n_J1_J2p = norm(J1-J2p,'fro')/ n_J1;

if n_J1_J2p > 1e-4; 
   warning('J1-J2p is too big');
end

function test4_npaa_3d
% Fine model
imdl = mk_common_model('n3r2');
f1mdl = imdl.fwd_model;
f1mdl.solve = @aa_fwd_solve;
f1mdl.jacobian = @aa_calc_jacobian;
f1mdl.system_mat = @aa_calc_system_mat;

% Create 2D FEM of all NODES with z=0
n2d = f1mdl.nodes( (f1mdl.nodes(:,3) == 0), 1:2);
e2d = delaunayn(n2d);
c_mdl = eidors_obj('fwd_model','2d','elems',e2d,'nodes',n2d);

% Create f_mdl with coarse2fine mapping
c2f= mk_coarse_fine_mapping( f1mdl, c_mdl );
f2mdl= f1mdl;
f2mdl.coarse2fine = c2f;

imdl.jacobian_bkgnd.value= .6;
img = calc_jacobian_bkgnd( imdl);

t=cputime;
J1= aa_calc_jacobian( f1mdl, img)*c2f;
fprintf('aa_calc - full: t=%f\n',  cputime-t  );
n_J1 = norm(J1,'fro');

t=cputime;
J2= aa_calc_jacobian( f2mdl, img);
fprintf('aa_calc - c2f: t=%f\n',  cputime-t  );

n_J1_J2 = norm(J1-J2,'fro')/ n_J1;
if n_J1_J2 > 1e-4; 
   warning('J1-J2 is too big');
end


t=cputime;
J2p= perturb_jacobian( f2mdl, img);
fprintf('perturb_j - c2f: t=%f\n',  cputime-t  );

n_J1_J2p = norm(J1-J2p,'fro')/ n_J1;
if n_J1_J2p > 1e-4; 
   warning('J1-J2p is too big');
end

