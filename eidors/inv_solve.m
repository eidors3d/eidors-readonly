function img = inv_solve( inv_model, data1, data2)
% INV_SOLVE: calculate imag from an inv_model and data
% 
% inv_solve can be called as
%     img= inv_solve( inv_model, data1, data2)
%   if inv_model.reconst_type = 'difference'
% or
%     img= inv_solve( inv_model, data )
%   if inv_model.reconst_type = 'static'
%
%   if inv_model.reconst_to = 'nodes' then output
%      img.node_data has output data
%   else   reconst_to = 'elems' (DEFAULT) then output to
%      img.elem_data has output data
%
% in each case it will call the inv_model.solve
%
% data      is a measurement data structure
% inv_model is a inv_model structure
% img       is an image structure
%           or a vector of images is data1 or data2 are vectors
%
% For difference EIT:
% data1      => difference data at earlier time (ie homogeneous)
% data2      => difference data at later time   (ie inhomogeneous)
%
% data can be:
%   - an EIDORS data object
%
%   - an M x S matrix, where M is the total number
%         of measurements expected by inv_model
%
%   - an M x S matrix, where M is n_elec^2
%        when not using data from current injection
%        electrodes, it is common to be given a full
%        measurement set.  For example, 16 electrodes give
%        208 measures, but 256 measure sets are common.
%        Data will be selected based on fwd_model.meas_select.
%
% If S > 1 for both data1 and data2 then the matrix sizes must be equal
%
% Parameters:
%   inv_model.inv_solve.select_parameters: indices of parameters to return
%                         DEFAULT: return all paramteres
%  Scale solution (to correct for amplitude or other defects)
%   inv_model.inv_solve.scale_solution.offset
%   inv_model.inv_solve.scale_solution.scale

% (C) 2005-2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

inv_model= eidors_model_params( inv_model );
opts = parse_parameters( inv_model );

eidors_msg('inv_solve: %s', inv_model.solve,2);

if opts.abs_solve
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end

   imgc= feval( inv_model.solve, inv_model, ...
               filt_data(inv_model,data1) );
else
   if nargin~=3;
      error('two data sets are required for a difference reconstruction');
   end

   % expand data sets if one is provided that is longer
   data_width= max(num_frames(data1), num_frames(data2));

   fdata1 = filt_data( inv_model, data1, data_width );
   fdata2 = filt_data( inv_model, data2, data_width );

   % TODO: Check if solver can handle being called with multiple data
   imgc= feval( inv_model.solve, inv_model, fdata1, fdata2);
end

img = eidors_obj('image', imgc );
% If we reconstruct with a different 'rec_model' then
%  put this into the img
if isfield(inv_model,'rec_model')
   img.fwd_model= inv_model.rec_model;
end


if opts.select_parameters;
   img.elem_data = img.elem_data( opts.select_parameters, :);
end;

if ~opts.reconst_to_elems
  img.node_data= img.elem_data;
  img = rmfield(img,'elem_data');
end

% Scale if required
try; img.elem_data = opts.offset + opts.scale * img.elem_data; end
try; img.node_data = opts.offset + opts.scale * img.node_data; end

function opts = parse_parameters( imdl );
   if  strcmp(imdl.reconst_type,'static') || ...
       strcmp(imdl.reconst_type,'absolute')
      opts.abs_solve = 1;
   elseif strcmp(imdl.reconst_type,'difference')
      opts.abs_solve = 0;
   else
      error('inv_model.reconst_type (%s) not understood', imdl.reconst_type); 
   end

   opts.select_parameters = [];
   try
      opts.select_parameters = imdl.inv_solve.select_parameters;
   end;

   opts.reconst_to_elems = 1;
   try; if strcmp( imdl.reconst_to, 'nodes' )
      opts.reconst_to_elems = 0;
   end; end
   
   opts.scale  = 1;
   try; opts.scale = imdl.inv_solve.scale_solution.scale; end

   opts.offset = 0;
   try; opts.offset = imdl.inv_solve.scale_solution.offset; end
 

% TODO: this code really needs to be cleaned, but not before eidors 3.4
function nf= num_frames(d0)
   if isnumeric( d0 )
      nf= size(d0,2);
   elseif d0(1).type == 'data';
      nf= size( horzcat( d0(:).meas ), 2);
   else
      error('Problem calculating number of frames. Expecting numeric or data object');
   end
   
% test for existance of meas_select and filter data
function d2= filt_data(inv_model, d0, data_width )
   if ~isnumeric( d0 )
       % we probably have a 'data' object

       d1 = [];
       for i=1:length(d0)
          if strcmp( d0(i).type, 'data' )
              d1 = [d1, d0(i).meas];
          else
              error('expecting an object of type data');
          end
       end

   else
      % we have a matrix of data. Hope for the best
      d1 = d0;
   end

   d1= double(d1); % ensure we can do math on our object

   if isfield(inv_model.fwd_model,'meas_select') && ...
     ~isempty(inv_model.fwd_model.meas_select)
      % we have a meas_select parameter that isn []

      meas_select= inv_model.fwd_model.meas_select;
      if     size(d1,1) == length(meas_select)
         d2= d1(meas_select,:);
      elseif size(d1,1) == sum(meas_select==1)
         d2= d1;
      else
         error('inconsistent difference data: (%d ~= %d). Maybe check fwd_model.meas_select',  ...
               d2_width, data_width);
      end
   else
      d2= d1;
   end

   if nargin==3 % expand to data width
      d2_width= size(d2,2);
      if d2_width == data_width
         % ok
      elseif d2_width == 1
         d2= d2(:,ones(1,data_width));
      else
         error('inconsistent difference data: (%d ~= %d)',  ...
               d2_width, data_width);
      end
   end

% Test code
function do_unit_test
   k=0; N=5; nd = 5;

   imdl = mk_common_model('d2c2',16);
   imdl = select_imdl( imdl, {'Choose NF=2.5'});
   mvx = linspace(-0.8,0.8,nd);
   [vh,vi] = simulate_movement(mk_image(imdl), [mvx;0*mvx;0.05+0*mvx]);

   img= inv_solve(imdl,vh,vi); img.show_slices.img_cols = 5;
   k=k+1; subplot(N,1,k); show_slices(img);
   unit_test_cmp('inv_solve: 1a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 1b',  std(img.elem_data), 1.5e-2, 1e-2);

   vhm= eidors_obj('data','nn','meas',vh);
   img= inv_solve(imdl,vhm,vi);
   unit_test_cmp('inv_solve: 2a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 2b',  std(img.elem_data), 1.5e-2, 1e-2);

   img= inv_solve(imdl,vh*ones(1,nd),vi);
   unit_test_cmp('inv_solve: 3a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 3b',  std(img.elem_data), 1.5e-2, 1e-2);

   vim= eidors_obj('data','nn','meas',vi);
   img= inv_solve(imdl,vhm,vim);
   unit_test_cmp('inv_solve: 4a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 4b',  std(img.elem_data), 1.5e-2, 1e-2);

   vhm= eidors_obj('data','nn','meas',vh*ones(1,nd));
   img= inv_solve(imdl,vhm,vi);
   unit_test_cmp('inv_solve: 5a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 5b',  std(img.elem_data), 1.5e-2, 1e-2);

   vhm(1)= eidors_obj('data','nn','meas',vh*ones(1,2));
   vhm(2)= eidors_obj('data','nn','meas',vh*ones(1,nd-2));
   img= inv_solve(imdl,vhm,vi);
   unit_test_cmp('inv_solve: 6a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 6b',  std(img.elem_data), 1.5e-2, 1e-2);

   vim(1)= eidors_obj('data','nn','meas',vi(:,1:3));
   vim(2)= eidors_obj('data','nn','meas',vi(:,4:end));
   img= inv_solve(imdl,vhm,vim);
   unit_test_cmp('inv_solve: 7a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 7b',  std(img.elem_data), 1.5e-2, 1e-2);

   im2 = imdl; im2.inv_solve.select_parameters = 1:5;
   img= inv_solve(im2,vh,vi);
   unit_test_cmp('inv_solve: 10', size(img.elem_data), [5,nd]);


   im2 = imdl;
   im2.inv_solve.scale_solution.offset = 1;
   im2.inv_solve.scale_solution.scale = 2;
   img= inv_solve(im2,vh,vi); img.show_slices.img_cols = 5;
   unit_test_cmp('inv_solve: 20a', mean(img.elem_data), 1.006, 2e-3);
   unit_test_cmp('inv_solve: 20b',  std(img.elem_data), 3e-2, 1e-2);
   
   im2.inv_solve.scale_solution.offset = 0;
   d = interp_mesh( imdl.fwd_model); d= sqrt(sum(d.^2,2));
   im2.inv_solve.scale_solution.scale = spdiags(1-d,0,length(d),length(d));
   img= inv_solve(imdl,vh,vi); img.show_slices.img_cols = 5;
   k=k+1; subplot(N,1,k); show_slices(img);

   im2 = select_imdl(imdl, {'Nodal GN dif'} );
   img= inv_solve(im2,vh,vi);
   unit_test_cmp('inv_solve: 30a', mean(img.node_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 30b',  std(img.node_data), 1.5e-2, 1e-2);

   im2.inv_solve.scale_solution.offset = 1;
   im2.inv_solve.scale_solution.scale = 2;
   img= inv_solve(im2,vh,vi); img.show_slices.img_cols = 5;
   unit_test_cmp('inv_solve: 31a', mean(img.node_data), 1.006, 2e-3);
   unit_test_cmp('inv_solve: 31b',  std(img.node_data), 3e-2, 1.5e-2);

   im2 = select_imdl(imdl, {'Basic GN abs'} );
   img= inv_solve(im2,vi(:,1));
   unit_test_cmp('inv_solve: 40a', mean(img.elem_data), 1.004, 1e-3);
   unit_test_cmp('inv_solve: 40b',  std(img.elem_data), 1.5e-2, 1e-2);


