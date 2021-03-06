function J= jacobian_adjoint( fwd_model, img)
% JACOBIAN_ADJOINT: J= jacobian_adjoint( img ) 
% Calculate Jacobian Matrix for current stimulation EIT
% J         = Jacobian matrix
% img.fwd_model = forward model
% img.elem_data = background for jacobian calculations
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'), do_unit_test, return, end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_ADJOINT with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end

img = data_mapper(img);
measurement = input_check(img);

org_params = img.current_params;
% all calcs use conductivity
img = convert_img_units(img, 'conductivity');

img.elem_data = check_elem_data(img);

fwd_model= img.fwd_model;

pp= fwd_model_parameters( fwd_model );

 J= do_jacobian_calculation( img, pp, fwd_model, org_params);
%J= do_jacobian_calculation_qr( img, pp, fwd_model);


% if ~strcmp(org_params,'conductivity')
%     J = apply_chain_rule(J, img, org_params);
%     if isfield(fwd_model, 'coarse2fine') && ...
%           size(img.elem_data,1)==size(fwd_model.coarse2fine,1)
%             J=J*fwd_model.coarse2fine;
%             nparam = size(fwd_model.coarse2fine,2);
%     end
% end

%restore img to original condition
if ~strcmp(org_params,'conductivity') || isfield(img, org_params)
    img = rmfield(img,'elem_data');
    img.current_params = [];
end

if ~strcmp(measurement, 'voltage')
   J = convert_measurement(J, img, measurement);
end

% calculate normalized Jacobian
if pp.normalize
   data= fwd_solve( img );
%  J= J ./ (data.meas(:)*ones(1,nparam));
   J= bsxfun(@rdivide, J, data.meas(:));
end

function J= do_jacobian_calculation( img, pp, fwd_model, org_params);

   s_mat= calc_system_mat( img );

   idx= 1:size(s_mat.E,1);
   idx( fwd_model.gnd_node ) = [];

   sv= zeros(pp.n_node, pp.n_stim );
   sv( idx,:) = left_divide(s_mat.E(idx,idx) , pp.QQ( idx,: ));

   zi2E= zeros(pp.n_elec, pp.n_node);
   % the minus below used to be missing
   zi2E(:, idx)= -pp.N2E(:,idx)/ s_mat.E(idx,idx) ;

   FC= system_mat_fields( fwd_model );
   do_c2f = false;
   if isfield(fwd_model,'coarse2fine') 
     do_c2f = true;
   end

   if do_c2f
      c2f = fwd_model.coarse2fine;
      if ~strcmp(org_params, 'conductivity')
         c2f = spdiag(chain_rule_factor(img, org_params)) * c2f;
      end
       
      DE = jacobian_calc(pp, zi2E, FC, sv, c2f);
      nparam= size(fwd_model.coarse2fine,2);
   else
      DE = jacobian_calc(pp, zi2E, FC, sv);
      nparam= pp.n_elem;
   end


   J = assemble_J(pp, nparam, DE, fwd_model.stimulation);
   if DEBUG
      Jo= assemble_J_old(pp, nparam, DE, fwd_model.stimulation);
      if norm(J-Jo,'fro')>1e-13;
         error('EIDORS:InternalValidation','assemble_J versions not match');
      end
   end
   
   % apply chain rule (already done for c2f)
   if ~do_c2f && ~strcmp(org_params, 'conductivity')
     J = apply_chain_rule(J, img, org_params);
   end

% New function which tries to qr decompose so we don't do
% extra calculations of extra stim patterns
function J= do_jacobian_calculation_qr( img, pp, fwd_model);
   s_mat= calc_system_mat( img );
   idx= 1:size(s_mat.E,1);
   idx( fwd_model.gnd_node ) = [];

   zi2E= zeros(pp.n_elec, pp.n_node);
   % the minus below used to be missing
   zi2E(:, idx)= -pp.N2E(:,idx)/ s_mat.E(idx,idx) ;

   FC= system_mat_fields( fwd_model );


   [Q,R] = qr(pp.QQ(idx,:),0);
   rnotzeros = any(R~=0,2);
   Q= Q(:,rnotzeros);
   R= R(rnotzeros,:);

if 0
   sv= zeros(pp.n_node, pp.n_stim );
   sv( idx,:) = left_divide(s_mat.E(idx,idx) , pp.QQ( idx,: ));
   jacobian_calc(pp, zi2E, FC, sv);
end

   sv= zeros(pp.n_node, sum(rnotzeros) );
   sv( idx,:) = s_mat.E(idx,idx) \ Q;
%t= cputime;
   DE = jacobian_calc(pp, zi2E, FC, sv);
%disp(cputime-t);
   DEr=  reshape(DE,size(sv,2),[]);
   Rt = R.';
%  pp.n_stim = size(sv,2);
%  Jn = assemble_J_old(pp, nparam, DE, fwd_model.stimulation);
   nparam= pp.n_elem;
   J = zeros( pp.n_meas, nparam );
   if 0
      idx=0;
      for j= 1:pp.n_stim
         meas_pat= fwd_model.stimulation(j).meas_pattern;
         nonzero = any(meas_pat,1);
         n_meas  = size(meas_pat,1);
   %     DEj = reshape( DE(:,j,:), pp.n_elec, [] );
   %     DEj = 0; for k = 1:size(R,1); DEj = DEj + full(R(k,j))* DE(k,:,:); end
   %     DEj = Rt(j,:) * DEr;
%        DEj = reshape( DEj, pp.n_elec, [] );
%        J( idx+(1:n_meas),: ) = meas_pat(:,nonzero)*DEj(nonzero,:);
         DEj = Rt(j,:) * reshape(DE(:,nonzero,:),size(sv,2),[]);
         DEj = reshape( DEj, sum(nonzero), [] );
         J( idx+(1:n_meas),: ) = meas_pat(:,nonzero)*DEj;
         idx= idx+ n_meas;
      end
   else
      stim = fwd_model.stimulation;
      J_out = arrayfun(@(s,j) multiplier(s,j,Rt,DE, pp), ...
           stim, 1:pp.n_stim, 'UniformOutput', false);
%             s.meas_pattern * reshape(Rt(j,:) * DEr, pp.n_elec,[]), ...
      J = vertcat(J_out{:});
   end

function Jj = multiplier(s,j,Rt,DE, pp);
%     Jj = s.meas_pattern * reshape(Rt(j,:) * DEr, pp.n_elec,[]);
      meas_pat = s.meas_pattern;
      nonzero = any(meas_pat,1);
      DEj = Rt(j,:) * reshape(DE(:,nonzero,:),size(DE,1),[]);
      DEj = reshape( DEj, sum(nonzero), [] );
      Jj = meas_pat(:,nonzero)*DEj;

function J = assemble_J_old(pp, nparam, DE, stim)
   J = zeros( pp.n_meas, nparam );
   idx=0;
   for j= 1:pp.n_stim
      meas_pat= stim(j).meas_pattern;
      n_meas  = size(meas_pat,1);
      DEj = reshape( DE(j,:,:), pp.n_elec, [] );
      J( idx+(1:n_meas),: ) = meas_pat*DEj;
      idx= idx+ n_meas;
   end

% Faster implementation of previous code
function J = assemble_J(pp, nparam, DE, stim)
   J_out = arrayfun(@(s,j) ...
           s.meas_pattern * reshape(DE(j,:,:),pp.n_elec,[]), ...
        stim, 1:pp.n_stim, 'UniformOutput', false);
   J = vertcat(J_out{:});

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
function DE = jacobian_calc(pp, zi2E, FC, sv, c2f)
d= pp.n_dims+1;
dfact= (d-1)*(d-2); % Valid for d<=3

do_c2f = ( nargin==5 );

zi2E_FCt = zi2E * FC';
FC_sv   = FC * sv;
n_stim = size(sv,2);

if ~do_c2f
   switch 2
   case 1
      DE= zeros(n_stim, pp.n_elec, pp.n_elem);
      for k= 1:pp.n_elem
          idx= (d-1)*(k-1)+1 : (d-1)*k;
          dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
%         DE(:,:,k)= dq.';
          dq = dq.';
          DE(:,:,k)= dq;
      end
   case 2
      DE= zeros(pp.n_elec, n_stim, pp.n_elem);
      for k= 1:pp.n_elem
          idx= (d-1)*(k-1)+1 : (d-1)*k;
          dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
          DE(:,:,k)= dq;
      end
      DE = permute(DE,[2,1,3]);
   case 3
      DE= zeros(pp.n_elem, pp.n_elec, n_stim);
      for k= 1:pp.n_elem
          idx= (d-1)*(k-1)+1 : (d-1)*k;
          dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
          DE(k,:,:)= dq;
      end
      DE = permute(DE,[3,2,1]);
   end
else
   DE= zeros(n_stim, pp.n_elec, size(c2f,2) );
   if 0 % Code is slower
      de= pp.n_elem * (d-1);
      for k= 1:size(c2f,2);
          chg_col = kron( c2f(:,k), ones(d-1,1));
          dDD_dEj = spdiags(chg_col,0, de, de);
          dq= zi2E_FCt * dDD_dEj * FC_sv;
          DE(:,:,k)= dq.';
      end
   else
      de= pp.n_elem * (d-1);
      for k= 1:size(c2f,2);
          ff = find( c2f(:,k) );
          lff= length(ff)*(d-1);
          ff1= ones(d-1,1) * ff(:)';
          ffd= (d-1)*ff1 + (-(d-2):0)'*ones(1,length(ff));
          dDD_dEj = spdiags(c2f(ff1,k), 0, lff, lff);
          dq= zi2E_FCt(:,ffd) * dDD_dEj * FC_sv(ffd,:);
          DE(:,:,k)= dq.';
      end
   end
end

function d = DEBUG; d=eidors_debug('query','jacobian_adjoint');

function dCond_dPhys = chain_rule_factor(img,org_params)
switch(org_params)
    case 'resistivity'
        dCond_dPhys = -img.elem_data.^2;
    case 'log_resistivity'
        dCond_dPhys = -img.elem_data;
    case 'log_conductivity'
        dCond_dPhys = img.elem_data;
    otherwise
        error('not implemented yet')
end


function J = apply_chain_rule(J, img, org_params)

dCond_dPhys = chain_rule_factor(img, org_params);
J = J*spdiag(dCond_dPhys);

function J = convert_measurement(J, img, measurement)
switch measurement
   case 'abs_voltage'
      img.fwd_model.measured_quantity = 'voltage';
      vv = fwd_solve(img);
      flip = sign(vv.meas);
      fctr = spdiags(flip, 0, length(flip), length(flip));    
   case 'log_voltage'
      img.fwd_model.measured_quantity = 'voltage';
      vv = fwd_solve(img);
      fctr = spdiag(1./vv.meas);
   case 'log10_voltage'
      img.fwd_model.measured_quantity = 'voltage';
      vv = fwd_solve(img);
      fctr = spdiag(1./(log(10)*vv.meas));
   case 'apparent_resistivity'
      fctr = apparent_resistivity_factor(img.fwd_model);
   case 'log_apparent_resistivity'
      img.fwd_model.measured_quantity = 'voltage';
      vv = fwd_solve(img);
      fctr = spdiag(log( diag(apparent_resistivity_factor(img.fwd_model))) ./ vv.meas);
   case 'log10_apparent_resistivity'
      img.fwd_model.measured_quantity = 'voltage';
      vv = fwd_solve(img);
      fctr = spdiag(log( diag(apparent_resistivity_factor(img.fwd_model))) ./(log(10)*vv.meas));
end
J = fctr*J;

function measurement = input_check(img)
if ~ismember(img.current_params, supported_params)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'JACOBIAN_ADJOINT',img.current_params);
end

try 
    measurement = img.fwd_model.measured_quantity;
catch
    measurement = 'voltage';
end
if ~ismember(measurement, supported_measurement)
    error('EIDORS:MeasurementNotSupported', '%s does not support %s',...
    'FWD_SOLVE_1ST_ORDER',measurement);
end

function elem_data = check_elem_data(img)
   elem_data = img.elem_data; 
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('jacobian_adjoin: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(img.fwd_model, 'coarse2fine');
     c2f = img.fwd_model.coarse2fine;
     sz_c2f = size(c2f);
     switch sz_elem_data(1)
       case sz_c2f(1); % Ok     
       case sz_c2f(2); elem_data = c2f * elem_data;
       otherwise; error(['jacobian_adjoint: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(img.fwd_model)
       error(['jacobian_adjoint: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end



function str = supported_params
    str = {'conductivity'
           'resistivity'
           'log_conductivity'
           'log_resistivity'};
        
function list = supported_measurement
   list = {'voltage'
           'log_voltage'
           'log10_voltage'
           'apparent_resistivity'
           'log_apparent_resistivity'
           'log10_apparent_resistivity'
           'abs_voltage'
           };
      
function do_unit_test
   J = jacobian_adjoint(mk_image(mk_common_model('a2c2',16),1));
   tst = [-0.000866389089603  -0.004906480672346   0.001753326247967];
   unit_test_cmp('Basic J', tst, J([2,22,222]), 1e-13);


   fmdl.type = 'fwd_model';
   fmdl.nodes = [0 0; 0 1; 1 1; 1 0; .5 .5];
   fmdl.elems = [1 2 5; 2 5 3; 3 5 4; 4 5 1];
   for i = 1:4
      fmdl.electrode(i).nodes = i;
      fmdl.electrode(i).z_contact = 1;
   end
   fmdl.stimulation.stim_pattern = sparse([1 0 -1 0]');
   fmdl.stimulation.meas_pattern = [1,-1, 0, 0
                                    0, 1,-1, 0
                                    0, 0, 1,-1
                                   -1, 0, 0, 1];
   fmdl.solve = 'eidors_default';
   fmdl.system_mat = 'eidors_default';
   fmdl.jacobian = @jacobian_adjoint;
   fmdl.normalize_measurements = 0;
   fmdl.gnd_node = 3;

% simple image   
%    show_fem(fmdl);
   J= jacobian_adjoint( mk_image( fmdl,1) );
   tst = [-1,0,0,-1;0,-1,-1,0;0,1,1,0;1,0,0,1];
   unit_test_cmp('Basic J', tst/2, J, 1e-13);
 
   fctr= apparent_resistivity_factor(fmdl);  
   elem_data_init= rand(size(fmdl.elems,1),1);
   meas = supported_measurement;
   phys = supported_params;
   j= 1; i = 1;
   img_test = mk_image(fmdl,1,phys{j});
   img_test.(phys{j}).elem_data = elem_data_init;
   img_test.fwd_model.measured_quantity = meas{i};
   J_init = jacobian_adjoint(img_test);
   v= fwd_solve(img_test);
   Voltage_estimated= v.meas;
   
   elem_data= [elem_data_init,1./elem_data_init,log(elem_data_init),-log(elem_data_init)];

   J_test(:,:,1,1)= J_init;
   J_test(:,:,1,2)= spdiag(1./(Voltage_estimated))*J_init;
   J_test(:,:,1,3)= spdiag(1./(log(10)*Voltage_estimated))*J_init;
   J_test(:,:,1,4)= fctr*J_init;
   J_test(:,:,1,5)= spdiag(log(diag(fctr))./Voltage_estimated)*J_init;
   J_test(:,:,1,6)= spdiag(log(diag(fctr))./(log(10)*Voltage_estimated))*J_init;
   J_test(:,:,1,7)= spdiag(sign(Voltage_estimated))*J_init;

   factor2= repmat(-elem_data_init'.^2,size(J_init,1),1);
   factor3= repmat(elem_data_init',size(J_init,1),1);
   factor4= repmat(-elem_data_init',size(J_init,1),1);
   for i=1:size(J_test,4);
      J_test(:,:,2,i)= factor2.*J_test(:,:,1,i);
      J_test(:,:,3,i)= factor3.*J_test(:,:,1,i);
      J_test(:,:,4,i)= factor4.*J_test(:,:,1,i);
   end

   tol= 1e-5;
   for j = 1:length(phys)
      img = mk_image(fmdl,1,phys{j});
      img.(phys{j}).elem_data = elem_data(:,j);
      for i = 1:length(meas)
         img.fwd_model.measured_quantity = meas{i};
         J = jacobian_adjoint(img);
         test= [ meas{i} ' vs ' phys{j} ];
         unit_test_cmp(test,J,J_test(:,:,j,i),tol);
      end
   end
