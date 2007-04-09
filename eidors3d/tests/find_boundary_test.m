function find_boundary_test
ok=1;

%2D Test:  
mdl2= mk_common_model('c2c',16);
bdy2= find_boundary(mdl2.fwd_model.elems);
bdy2= sort_boundary(bdy2);
bdyc= sort_boundary(mdl2.fwd_model.boundary);

ok= match(bdy2,bdyc,ok,'2D test');



function bdy= sort_boundary(bdy)
   bdy = sort(bdy,2);
   [jnk,idx] =sort(bdy(:,1));
   bdy = bdy(idx,:);

function ok= match( pat1, pat2, ok, descr)
   if ~all(pat1(:) == pat2(:))
      ok=0;
keyboard
      eidors_msg('find_bounday_test: fail %s',descr,1);
   end
