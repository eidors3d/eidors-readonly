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

% find which nodes are in contact
function pp= base_n_contact(img,pp);
   fmdl = img.fwd_model;
   base= fmdl.remesh_solver.base_nodes(:);
   elems = fmdl.elems;
   gone = ismember(elems, base);
   pp.base_elems = all(gone,2);
   % Eliminate elems only connected to base
   elems(pp.base_elems,:) = [];
   pp.contact_nodes = intersect(elems(:), base);
   pp.base_nodes = setxor(base,pp.contact_nodes);
   

% base nodes (with ground), without, contact
function pp = segment_model(img);
   fmdl = img.fwd_model;
   params = fmdl.remesh_solver;
   pp= fwd_model_parameters(fmdl);

   Nnode= num_nodes(img);
   pp= base_n_contact(img,pp);

   base_elecs= params.base_elecs(:);
   extr_elecs= params.extr_elecs(:);
   pp.baseg= [pp.base_nodes; Nnode + base_elecs];

   pp.base = pp.baseg;
   pp.base(pp.base == fmdl.gnd_node) = [];

   pp.cont = pp.contact_nodes(:);

   pp.extr= [(1:num_nodes(img))';
              Nnode + extr_elecs];
   pp.extr(pp.base_nodes) = [];

   [~,pp.idx,~] = intersect(pp.extr,pp.cont);

% TODO: There are more efficient ways of
%  caching base model.
function [vA,Ab,bAb,B] = base_block(img,pp);
   nbe = ~pp.base_elems;
   img.elem_data(nbe,:) = [];  
   img.fwd_model.elems(nbe,:) = [];  
   E = calc_system_mat(img); E=E.E;

   A = E(pp.base,pp.base);
   B = E(pp.base,pp.cont);
   iA= pp.QQ(pp.base,:);
   [vA,Ab,bAb] = eidors_cache( ...
       @base_block_calc, {A,B,iA});



function [vA,Ab,bAb]= base_block_calc(A,B,iA);
%  vA = left_divide(A,iA);
%  Ab = left_divide(A,B);
% Solve once so LU only needs to be done once
   vA_Ab = left_divide(A,[iA,B]);
   LiA = size(iA,2);
   vA = vA_Ab(:,1:LiA);
   Ab = vA_Ab(:,LiA+1:end);
   bAb= B'*Ab;


% TODO: Calcs full matrix
function D = inner_sys_mat_old(img,pp);
   E = calc_system_mat(img); E=E.E;
   D = E(pp.extr,pp.extr);


% TODO: Avoid calculating full system mat
function D = inner_sys_mat(img,pp);
   gone= all( ~ismember( ...
        img.fwd_model.elems,pp.extr),2);
   img.elem_data(gone,:) = [];
   img.fwd_model.elems(gone,:) = [];
   E = calc_system_mat(img); E=E.E;
   D = E(pp.extr,pp.extr);


function [vD,Db]= extr_block(img,bAb,B,pp);
   DmA = inner_sys_mat(img,pp);
   ii = pp.idx;

   iD= pp.QQ(pp.extr,:);
   DmA(ii,ii)= DmA(ii,ii) - bAb;

   switch 2;
      case 1;
   vD = left_divide(DmA,iD);

   Bx= sparse(size(vD,1),size(B,1));
   Bx(ii,:) = B';
   Db = left_divide(DmA,Bx);
      case 2;
   iD_Bx = iD;
   LiD = size(iD,2);
   LBx = size(B,1);
   iD_Bx(ii,LiD+(1:LBx)) = B';
   vD_Db = left_divide(DmA, iD_Bx);
   vD = vD_Db(:,1:LiD);
   Db = vD_Db(:,LiD+1:end);
      case 3;
   DmAi = inv(DmA);      
   vD = DmAi * iD;
   Db = DmAi(:,ii)*B';
      end
   

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
img.fwd_model.remesh_solver.base_nodes = 5:num_nodes(fmdl);
img.fwd_model.remesh_solver.base_elecs = 1:4;
img.fwd_model.remesh_solver.extr_elecs = 5;
vn=remesh_solver(img);
vf=fwd_solve(img);
subplot(122);plot(vf.volt - vn.volt); box off;
unit_test_cmp('a2C2',vf.volt,vn.volt,1e-12);
