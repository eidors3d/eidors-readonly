function elemcur = calc_elem_current( img, vv )
% calc_elem_current: calculate current vector in each FEM element
%
% e_curr = calc_elem_current( img, vv )
%   img -> img object 
%   volt-> voltage on nodes if not specified, img is solved
%      via fwd_solve
%   e_curr = current in each element [N_elems x N_dims]


% (C) 2017 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

dims = size(img.fwd_model.nodes,2);

if nargin==1; % We need to calculate
   if isfield(img,'elem_data')
      img.fwd_solve.get_all_meas = 1;
      vh = fwd_solve(img);
      vv = vh.volt(:,1);
   elseif isfield(img,'node_data');
      vv = img.node_data(:,1);
      error('show_current: cannot interpolate conductivity onto elements (yet)');
   else
      error('show_current: one parameter provided, and cannot solve for voltages');
   end
end 

Nel = size(img.fwd_model.elems,1);
elemcur = zeros(Nel,dims);
% Calc field as I = sigma E
%V1 = V0 + Ex*x1 + Ey*y1   [ 1 x1 y1 ] [ V0 ]
%V2 = V0 + Ex*x2 + Ey*y2 = [ 1 x2 y2 ]*[ Ex ]
%V3 = V0 + Ex*x3 + Ey*y    [ 1 x3 y3 ] [ Ey ]
oo = ones(dims+1,1);
for i=1:Nel
  idx = img.fwd_model.elems(i,:);
  nod = img.fwd_model.nodes(idx,:);
  if dims ==2
     VE  = ([oo, nod])\fix_dim(vv(idx));
  else
     VE  = ([oo, nod])\vv(idx);
  end
  elemcur(i,:) = img.elem_data(i,1)*VE(2:end)';
%  elemcur(i+1,:) = (reshape(img.elem_data(i,1,:,:),dims,dims)*VE(2:end))';
end

% In case it is the wrong vector shape
function vv = fix_dim(vv)
    if size(vv,1) == 1
        vv = vv';
    end

function do_unit_test
   fmdl.nodes = [0,0;0,1;1,0;1,1];
   fmdl.elems = [1,2,3;2,3,4];
   fmdl.electrode(1).nodes = [1,2]; fmdl.electrode(1).z_contact = 0.01;
   fmdl.electrode(2).nodes = [3,4]; fmdl.electrode(2).z_contact = 0.01;
   fmdl.gnd_node = 1;
   fmdl.stimulation(1).stim_pattern = [1;-1];
   fmdl.stimulation(1).meas_pattern = [1,-1];
   fmdl.solve = @fwd_solve_1st_order;
   fmdl.system_mat = @system_mat_1st_order;
   fmdl.type = 'fwd_model';
   fmdl.normalize_measurements= 0;
   img = mk_image(fmdl,[1,1]); 
   img.fwd_solve.get_all_meas = 1;

   e_curr = calc_elem_current(img);
   unit_test_cmp('simple_mdl:', e_curr, [-1,0;-1,0],1e-12);


   imdl= mk_common_model('d2c2',8);
   img = calc_jacobian_bkgnd( imdl );
   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve(img);
   show_current(img, vh.volt(:,1));
   e_curr = calc_elem_current(img);
   unit_test_cmp('d2c2:', e_curr([1,10,100],:), ...
          [2.422593061890268, -0.920998260630422;
           2.887551382948032, -1.051869372020626;
           1.349507157073426, -0.842871247870084],1e-12);
