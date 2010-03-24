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
%   else   reconst_to = 'elems' then output to
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
% If S > 1 for both data1 and data2 then the values must be equal
%
% Parameters:
%   inv_model.inv_solve.select_parameters: indices of parameters to return
%                         DEFAULT: return all paramteres

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% COMMENT: There seems to be no general way to cache
%       inv_model parameters. Thus, each algorithm needs
%       to do it separately. It may be possible to cache
%       one-step inverse matrices, but that is not done. 

inv_model= eidors_model_params( inv_model );
eidors_msg('inv_solve: %s', inv_model.solve,2);


if     strcmp(inv_model.reconst_type,'static') || ...
       strcmp(inv_model.reconst_type,'absolute')
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end

   imgc= feval( inv_model.solve, inv_model, ...
               filt_data(inv_model,data1) );

elseif strcmp(inv_model.reconst_type,'difference')
   if nargin~=3;
      error('two data sets are required for a difference reconstruction');
   end

   % expand data sets if one is provided that is longer
   data_width= max(num_frames(data1), num_frames(data2));

   fdata1 = filt_data( inv_model, data1, data_width );
   fdata2 = filt_data( inv_model, data2, data_width );

   % Check if solver can handle being called with multiple data
   imgc= feval( inv_model.solve, inv_model, fdata1, fdata2);
else
   error('inv_model.reconst_type not understood'); 
end

img = eidors_obj('image', imgc );
% If we reconstruct with a different 'rec_model' then
%  put this into the img
if isfield(inv_model,'rec_model')
   img.fwd_model= inv_model.rec_model;
end

try; if length(inv_model.inv_solve.select_parameters)>0
   img.elem_data = img.elem_data( ...
               inv_model.inv_solve.select_parameters,:);
end; end

try % move elem_data to nodes if required
   if inv_model.reconst_to == 'nodes'
      img.node_data= img.elem_data;
      img = rmfield(img,'elem_data');
   end
end

function nf= num_frames(d0)
   if isnumeric( d0 )
      nf= size(d0,2);
   elseif d0(1).type == 'data';
      nf= length(d0);
   else
      error('Problem calculating number of frames');
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

   if isfield(inv_model.fwd_model,'meas_select');
      % we have a meas_select parameter

      meas_select= inv_model.fwd_model.meas_select;
      if     size(d1,1) == length(meas_select)
         d2= d1(meas_select,:);
      elseif size(d1,1) == sum(meas_select==1)
         d2= d1;
      else
         error('data size does not match meas_select');
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

