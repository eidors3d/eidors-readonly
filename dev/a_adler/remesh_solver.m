function data =remesh_solver(img)
% img.fwd_model 
%    .remesh_solver.base_nodes= list
%    .remesh_solver.contact_nodes = list
%    .remesh_solver.base_elecs= list
%       CEM electrodes attached to the base
%    .remesh_solver.extr_elecs= list
%       CEM electrodes attached to the extr

% (C)2021 Andy Adler. Licenced under GPL v2 or v3

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test;return;end

pp = segment_model(img);
[vA,Ab,bAb,B] = base_block(img,pp);
[vD,Db]       = extr_block(img,bAb,B,pp);

vI = vD - Db*vA;
vO = vA - Ab*vI(pp.idx,:);
v([pp.base;pp.extr],:) = [vO; vI];

data = pack_output(img,v);

% base nodes (with ground), without, contact
function pp = segment_model(img);
   fmdl = img.fwd_model;
   params = fmdl.remesh_solver;
   pp= fwd_model_parameters(fmdl);

   Nnode= num_nodes(img);
   base_nodes= params.base_nodes(:);
   base_elecs= params.base_elecs(:);
   extr_elecs= params.extr_elecs(:);
   pp.baseg= [base_nodes; Nnode + base_elecs];

   pp.base = pp.baseg;
   pp.base(pp.base == fmdl.gnd_node) = [];

   pp.cont = fmdl.remesh_solver.contact_nodes(:);

   pp.extr= [(1:num_nodes(img))';
              Nnode + extr_elecs];
   pp.extr(base_nodes) = [];

   [~,pp.idx,~] = intersect(pp.extr,pp.cont);

% TODO: Cache on base model
function [vA,Ab,bAb,B] = base_block(img,pp);

   E = calc_system_mat(img); E=E.E;
   A = E(pp.base,pp.base);
   B = E(pp.base,pp.cont);
   iA= pp.QQ(pp.base,:);

   vA = left_divide(A,iA);
   Ab = left_divide(A,B);
   bAb= B'*Ab;

function E = inner_sys_mat_old(img,pp);
   E = calc_system_mat(img); E=E.E;

function E = inner_sys_mat(img,pp);
   E = calc_system_mat(img); E=E.E;



% TODO: Avoid calculating full system mat
function [vD,Db]= extr_block(img,bAb,B,pp);
   E = inner_sys_mat(img,pp);

   iD= pp.QQ(pp.extr,:);
   DmA = E(pp.extr,pp.extr);
   DmA(pp.idx,pp.idx)= DmA(pp.idx,pp.idx) - bAb;
   vD = left_divide(DmA,iD);

   Bx= zeros(size(vD,1),size(B,1));
   Bx(pp.idx,:) = B';
   Db = left_divide(DmA,Bx);

function data = pack_output(img,v)
   fmdl = img.fwd_model;
   pp = fwd_model_parameters(fmdl);

   idx = find(any(pp.N2E));
   data.volt = v;
   data.time= NaN; % unknown
   data.name= 'solved by block solver';

   v_els= pp.N2E(:,idx) * v(idx,:);
   idx=0;
   stim = fmdl.stimulation;
   for i=1:length(stim)
      meas_pat= stim(i).meas_pattern;
      n_meas  = size(meas_pat,1);
      vv( idx+(1:n_meas),: ) = meas_pat*v_els(:,i);
      idx= idx+ n_meas;
   end
   data.meas= vv;


function do_unit_test
fmdl= mk_common_model('a2C2',4);
fmdl=fmdl.fwd_model;
fmdl.elems(1:4,:) = [];
fmdl.gnd_node = 15;
fmdl = remove_unused_nodes(fmdl);
fmdl.electrode(5) = struct('nodes',1:3,'z_contact',1e-3);
fmdl.boundary = find_boundary(fmdl);
fmdl.stimulation = stim_meas_list([1,2,3,4;1,3,3,5],5);
img = mk_image(fmdl,1);
img.fwd_solve.get_all_nodes = true;
subplot(121);show_fem(img,[0,0,2.012]);
img.fwd_model.remesh_solver.base_nodes = 13:num_nodes(fmdl);
img.fwd_model.remesh_solver.contact_nodes = 5:12;
img.fwd_model.remesh_solver.base_elecs = 1:4;
img.fwd_model.remesh_solver.extr_elecs = 5;
vn=remesh_solver(img);
vf=fwd_solve(img);
subplot(122);plot(vf.volt - vn.volt); box off;
unit_test_cmp('a2C2',vf.volt,vn.volt,1e-12);
