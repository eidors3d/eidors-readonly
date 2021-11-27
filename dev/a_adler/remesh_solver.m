function data =nov26_remesh_solver(img)
% img.fwd_model 
%    .remesh_solver.base_nodes= list
%    .remesh_solver.base_elecs= list
%       CEM electrodes attached to the base

% (C)2021 Andy Adler. Licenced under GPL v2 or v3

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test;return;end

[vA,Ab,bAb,B,base,extr] = base_block(img);
[vD,Db]       = extr_block(img,extr,bAb,B);

vI = vD - Db*vA;
vO = vA - Ab*vI;
v([base(:);extr(:)],:) = [vO; vI];

data = pack_output(img,v);


function [vA,Ab,bAb,B,base,extr] = base_block(img);
   fmdl = img.fwd_model;
   pp= fwd_model_parameters(fmdl);
   E = calc_system_mat(img);       E=E.E;

   Nnode= num_nodes(img);
   base = [fmdl.remesh_solver.base_nodes(:);
     Nnode+fmdl.remesh_solver.base_elecs(:)];
   extr = 1:size(E,1); extr(base) = [];
   base(base == fmdl.gnd_node) = [];

   A = E(base,base);
   B = E(base,extr);
   iA= pp.QQ(base,:);

   vA = left_divide(A,iA);
   Ab = left_divide(A,B);
   bAb= B'*Ab;

function [vD,Db] = extr_block(img,extr,bAb,B);
   pp= fwd_model_parameters(img.fwd_model);
   E = calc_system_mat(img);       E=E.E;

   iD= pp.QQ(extr,:);
   D = E(extr,extr);
   DmA= D - bAb;
   vD = left_divide(DmA,iD);
   Db = left_divide(DmA,B');

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
fmdl= mk_common_model('a2C2',4); fmdl=fmdl.fwd_model;
fmdl.elems(1:4,:) = [];
fmdl = remove_unused_nodes(fmdl);
fmdl.electrode(5) = struct('nodes',1:3,'z_contact',1e-3);
fmdl.boundary = find_boundary(fmdl);
fmdl.gnd_node = 15;
fmdl.stimulation = stim_meas_list([1,2,3,4;1,3,3,5],5);
img = mk_image(fmdl,1);
img.fwd_solve.get_all_nodes = true;
subplot(121);show_fem(img,[0,0,2.012]);
img.fwd_model.remesh_solver.base_nodes = 13:num_nodes(fmdl);
img.fwd_model.remesh_solver.base_elecs = 1:4;
vn=remesh_solver(img);
vf=fwd_solve(img);
subplot(122);plot(vf.volt - vn.volt); box off;
unit_test_cmp('a2C2',vf.volt,vn.volt,1e-12);
