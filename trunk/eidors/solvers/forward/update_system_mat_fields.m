function FC1 = update_system_mat_fields( fwd_model0, fwd_model1 )
% SYSTEM_MAT_FIELDS: fields (elem to nodes) fraction of system mat
% FC= system_mat_fields( fwd_model )
% input: 
%   fwd_model = forward model
% output:
%   FC:        s_mat= C' * S * conduct * C = FC' * conduct * FC;

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model0) && strcmp(fwd_model0,'UNIT_TEST'); do_unit_test; return; end

copt.cache_obj = mk_cache_obj(fwd_model0);
copt.fstr = 'system_mat_fields';
copt.log_level = 4;
FC0 = eidors_cache(@system_mat_fields,{fwd_model0},copt ); % this should be a cache hit

t = tic;
copt.cache_obj = mk_cache_obj(fwd_model1);
copt.fstr = 'update_system_mat_fields';
dFC = eidors_cache(@calc_update_system_mat_fields,{fwd_model0, fwd_model1},copt );

% updated system_mat_fields
FC1 = FC0 + dFC;
t = toc(t);

% now we fake out the caching system by telling it that we are doing
% system_mat_fields(fwd_model1) with the new model, calls to system_mat_fields
% will get the cached value from this function
fstr = 'system_mat_fields';
cache_obj = mk_cache_obj(fwd_model1);
eidors_obj('set-cache', cache_obj, fstr, {FC1}, t);


% only cache stuff which is really relevant here
function cache_obj = mk_cache_obj(fwd_model)
   cache_obj.elems       = fwd_model.elems;
   cache_obj.nodes       = fwd_model.nodes;
   try
   cache_obj.electrode   = fwd_model.electrode; % if we have it
   end
   cache_obj.type        = 'fwd_model';
   cache_obj.name        = ''; % it has to have one

function dFC= calc_update_system_mat_fields( fwd_model0, fwd_model1 )
   p0= fwd_model_parameters( fwd_model0, 'skip_VOLUME' );
   p1= fwd_model_parameters( fwd_model1, 'skip_VOLUME' );
   d0= p0.n_dims+0;
   d1= p0.n_dims+1;
   e= p0.n_elem;
   n= p0.n_node;

   % find moved nodes, then the elements affected by those moved nodes
   [dn,~] = find(abs(fwd_model0.nodes - fwd_model1.nodes) > eps);
   dn = unique(dn);
   [de,~] = find(ismember(fwd_model1.elems,dn)); % our modified nodes touched these elements
   assert(all(all(fwd_model0.elems == fwd_model1.elems)), 'expected fmdl0 and fmdl1 to have the same element structure');
   assert(all(de <= e), 'invalid elem# found for delta nodes');

   FFjidx= floor([0:d0*e-1]'/d0)*d1*ones(1,d1) + ones(d0*e,1)*(1:d1);
   FFiidx= [1:d0*e]'*ones(1,d1);
   FFdata= zeros(d0*e,d1);
   dfact = (d0-1)*d0;
   for j=de(:)'
     a0=  inv([ ones(d1,1), p0.NODE( :, p0.ELEM(:,j) )' ]);
     a1=  inv([ ones(d1,1), p1.NODE( :, p0.ELEM(:,j) )' ]);
     idx= d0*(j-1)+1 : d0*j;
     FFdata(idx,1:d1)= a1(2:d1,:)/ sqrt(dfact*abs(det(a1))) - a0(2:d1,:)/ sqrt(dfact*abs(det(a0)));
   end %for j=1:ELEMs

if 0 % Not complete electrode model
   FF= sparse(FFiidx,FFjidx,FFdata);
   CC= sparse((1:d1*e),p0.ELEM(:),ones(d1*e,1), d1*e, n);
else
   % TODO this could be better: we could look to see which, if any electrodes
   % have been modified and only update those ones... currently the
   % implementation here is pretty brain dead on the assumption this part is
   % fast
   % TODO currently can't handle electrode node number changes
   [F2data0,F2iidx0,F2jidx0, C2data0,C2iidx0,C2jidx0] = compl_elec_mdl(fwd_model0,p0,dn);
   [F2data1,F2iidx1,F2jidx1, C2data1,C2iidx1,C2jidx1] = compl_elec_mdl(fwd_model1,p1,dn);

   FF= sparse([FFiidx(:);  F2iidx0(:);  F2iidx1(:)],...
              [FFjidx(:);  F2jidx0(:);  F2jidx1(:)],...
              [FFdata(:); -F2data0(:); +F2data1(:)]);
   
   CC= sparse([(1:d1*e)';     C2iidx0(:)], ...
              [p0.ELEM(:);    C2jidx0(:)], ...
              [ones(d1*e,1);  C2data0(:)]);
end

dFC= FF*CC;

% Add parts for complete electrode model
function [FFdata,FFiidx,FFjidx, CCdata,CCiidx,CCjidx] = ...
             compl_elec_mdl(fwd_model,pp,dn)
   d0= pp.n_dims;
   FFdata= zeros(0,d0);
   FFd_block= sqrtm( ( ones(d0) + eye(d0) )/6/(d0-1) ); % 6 in 2D, 12 in 3D 
   FFiidx= zeros(0,d0);
   FFjidx= zeros(0,d0);
   FFi_block= ones(d0,1)*(1:d0);
   CCdata= zeros(0,d0);
   CCiidx= zeros(0,d0);
   CCjidx= zeros(0,d0);
  
   sidx= d0*pp.n_elem;
   cidx= (d0+1)*pp.n_elem;
   for i= 1:pp.n_elec
      eleci = fwd_model.electrode(i);
      % contact impedance zc is in [Ohm.m] for 2D or [Ohm.m^2] for 3D
      zc=  eleci.z_contact;
      if any(find(ismember(eleci.nodes,dn)))
         [bdy_idx, bdy_area] = find_electrode_bdy( ...
             pp.boundary, fwd_model.nodes, eleci.nodes );
             % bdy_area is in [m] for 2D or [m^2] for 3D
   
         for j= 1:length(bdy_idx);
            bdy_nds= pp.boundary(bdy_idx(j),:);
   
            % 3D: [m^2]/[Ohm.m^2] = [S]
            % 2D: [m]  /[Ohm.m]   = [S]
            FFdata= [FFdata; FFd_block * sqrt(bdy_area(j)/zc)];
            FFiidx= [FFiidx; FFi_block' + sidx];
            FFjidx= [FFjidx; FFi_block  + cidx];
   
            CCiidx= [CCiidx; FFi_block(1:2,:) + cidx];
            CCjidx= [CCjidx; bdy_nds ; (pp.n_node+i)*ones(1,d0)];
            CCdata= [CCdata; [1;-1]*ones(1,d0)];
            sidx = sidx + d0;
            cidx = cidx + d0;
         end
      else
         [bdy_idx] = find_electrode_bdy( ...
             pp.boundary, fwd_model.nodes, eleci.nodes );
             % bdy_area is in [m] for 2D or [m^2] for 3D
   
         for j= 1:length(bdy_idx);
            bdy_nds= pp.boundary(bdy_idx(j),:);
   
            % 3D: [m^2]/[Ohm.m^2] = [S]
            % 2D: [m]  /[Ohm.m]   = [S]
            FFdata= [FFdata; FFd_block * 0];
            FFiidx= [FFiidx; FFi_block' + sidx];
            FFjidx= [FFjidx; FFi_block  + cidx];
   
            CCiidx= [CCiidx; FFi_block(1:2,:) + cidx];
            CCjidx= [CCjidx; bdy_nds ; (pp.n_node+i)*ones(1,d0)];
            CCdata= [CCdata; [1;-1]*ones(1,d0)];
            sidx = sidx + d0;
            cidx = cidx + d0;
         end
      end
      
   end

function do_unit_test
   disp('running system_mat_fields UNIT_TEST');
   system_mat_fields('UNIT_TEST'); % we depend explicitly on system_mat_fields... check it works properly!

   eidors_cache('clear_name', 'system_mat_fields');
   eidors_cache('clear_name', 'update_system_mat_fields');
   disp('running update_system_mat_fields UNIT_TEST');
   imdl=  mk_geophysics_model('h2e',32);
   ndim = size(imdl.fwd_model.nodes,2);
   for i = 1:10
      fmdl0 = imdl.fwd_model;
      fmdl0.nodes(1,:) = fmdl0.nodes(1,:) + rand(1,ndim)*1e-10; % defeat cache
      t=tic; FC0 = system_mat_fields(fmdl0); t0(i)=toc(t);
   
      % perturb nodes
      fmdl1 = fmdl0;
      nn = fmdl1.electrode(1).nodes;
      nn = [nn(1); nn(end)];
      fmdl1.nodes(nn,:) = fmdl1.nodes(nn,:) + 1e-4+rand(2,ndim)*1e-10; % perturb
      t=tic; FC1a = update_system_mat_fields(fmdl0, fmdl1); t1(i)=toc(t);
      t=tic; FC1b = system_mat_fields(fmdl1); t2(i)=toc(t);
      eidors_cache('clear_name', 'system_mat_fields');
      t=tic; FC1c = system_mat_fields(fmdl1); t3(i)=toc(t);
   end
   fprintf('system_mat_fields(fmdl0) = %0.2f sec\n',mean(t0));
   fprintf('update_system_mat_fields(fmdl0,fmdl1) = %0.3f sec [faster?]\n',mean(t1));
   fprintf('system_mat_fields(fmdl1) = %0.3f sec [cache hit?]\n',mean(t2));
   fprintf('system_mat_fields(fmdl1+delta) = %0.3f sec [recalculate]\n',mean(t3));
   unit_test_cmp('delta FC is fast                  ',mean(t1) < mean(t0)/2,1);
   unit_test_cmp('cache trick is really fast        ',mean(t2./t0) < 0.015,1);
   unit_test_cmp('cache was cleared before recalc   ',mean(t3) > mean(t0)*0.9,1); % did actually clear cache
   unit_test_cmp('cache trick returns correct result',FC1a,FC1b);
   unit_test_cmp('delta FC gives same result        ',FC1a,FC1c,10*eps);
   fprintf('speed-up: %0.2f\n',mean(t0./t1));
   if(sum(sum((FC1a - FC1c).^2)) > 0)
      err = FC1a - FC1c
   end

   if 0
      disp('---- profiling ----');
      fmdl1.nodes(nn,:) = fmdl1.nodes(nn,:) + rand(2,ndim)*1e-3; % perturb
      t=tic; FC0 = system_mat_fields(fmdl0); t0=toc(t);
      profile clear;
      profile on;
      t = tic; FC1a = update_system_mat_fields(fmdl0, fmdl1); t1=toc(t);
      profile off;
      profview;
      profsave(profile('info'),'profile');
      fprintf('update_system_mat_fields(fmdl0,fmdl1) = %0.3f sec [profiled]\n',t1);
   end
