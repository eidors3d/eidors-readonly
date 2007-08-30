function find_boundary_test
ok=1;

%2D Test:  
mdl = mk_common_model('c2c',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

ok= match(bdy,bdyc,ok,'2D test');

%3D Test:  
mdl = mk_common_model('n3r2',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

ok= match(bdy,bdyc,ok,'3D test n3r2');

%3D Test:  
mdl = mk_common_model('a3cr',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

% Bug is expected with old code in mk_circ_model - now uses find_boundary
ok= match(bdy,bdyc,ok,'3D test a3c2 - expected bug');


function bdy= sort_boundary(bdy)
   bdy = sort(bdy,2);
   bdy = sortrows(bdy);

function ok= match( pat1, pat2, ok, descr)
   if ~all(pat1(:) == pat2(:))
      ok=0;
      eidors_msg('find_bounday_test: fail %s',descr,1);
   end
