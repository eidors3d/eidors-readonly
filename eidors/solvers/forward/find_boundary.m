function [srf, idx] = find_boundary(simp);
% [srf, idx] = find_boundary(simp);
%
%Caclulates the boundary faces of a given 3D volume.
%Usefull in electrode assignment.
%
%srf  =  array of elements on each boundary simplex
%        boundary simplices are of 1 lower dimention than simp
%idx  =  index of simplex to which each boundary belongs
%simp = The simplices matrix

% $Id$

if isstr(simp) && strcmp(simp,'UNIT_TEST'); do_unit_test; return; end
if isstruct(simp) && strcmp(simp.type,'fwd_model'); simp= simp.elems; end

wew = size(simp,2) - 1;

if wew==3 || wew==2
   [srf,idx]= find_2or3d_boundary(simp,wew);
else
   eidors_msg('find_boundary: WARNING: not 2D or 3D simplices',1);
   srf=[]; return;
end

% sort surfaces. If there is more than one, its not on the boundary
function [srf,idx]= find_2or3d_boundary(simp,dim);
   if size(simp,1) < 4e9 % max of uint32
      % convert to integer to make sort faster
      simp = uint32( simp );
   end
   localface = nchoosek(1:dim+1,dim);
   srf_local= simp(:,localface');
   srf_local= reshape( srf_local', dim, []); % D x 3E
   srf_local= sort(srf_local)'; % Sort each row
   [sort_srl,sort_idx] = sortrows( srf_local );

   % Fine the ones that are the same
   first_ones =  sort_srl(1:end-1,:);
   next_ones  =  sort_srl(2:end,:);
   same_srl = find( all( first_ones == next_ones, 2) );

   % Assume they're all different. then find the same ones
   diff_srl = logical(ones(size(srf_local,1),1));
   diff_srl(same_srl) = 0;
   diff_srl(same_srl+1) = 0;

   srf= sort_srl( diff_srl,: );
   idx= sort_idx( diff_srl);
   idx= ceil(idx/(dim+1));

function do_unit_test
ok=1;

%2D Test:  
mdl = mk_common_model('c2c',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

ok= match(bdy,bdyc,ok,'2D test');

%3D Test:  
mdl = mk_common_model('n3r2',[16,2]);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

ok= match(bdy,bdyc,ok,'3D test n3r2');

%3D Test:  
mdl = mk_common_model('a3cr',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

ok= match(bdy,bdyc,ok,'3D test a3c2');

%3D Test:  
mdl = mk_common_model('b3cr',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

ok= match(bdy,bdyc,ok,'3D test b3c2');

function bdy= sort_boundary(bdy)
   bdy = sort(bdy,2);
   bdy = sortrows(bdy);

function ok= match( pat1, pat2, ok, descr)
    ok =  all(pat1(:) == pat2(:));
    fprintf('find_bounday_test: ok=%d\n',ok);

% function ok= match( pat1, pat2, ok, descr)
%    if ~all(pat1(:) == pat2(:))
%       ok=0;
%       eidors_msg('find_bounday_test: fail %s',descr,1);
%    end


