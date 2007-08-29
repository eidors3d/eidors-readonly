function dva = calc_difference_data( data1, data2, fwd_model)
% CALC_DIFFERENCE_DATA: calculate difference data between
%   two eidors data vectors
% dva = calc_difference_data( meas1, meas2, fwd_model)
%   
%  dva   = Matrix n_meas x n_time_steps of difference meas
%  meas1 = measurement object (or matrix) at time1 (homogeneous)
%  meas2 = measurement object (or matrix) at time2 (inhomogeneous)
%  fwd_model (optional, if provided in meas1 and meas2)
%
% This code appears simple, but there are a number of tricks
%  to remember, so it is best to factor it out. Issues are
%  1) normalize_data, 2) remove zero meas from adjacent systems,
%  3) allow both raw data and eidors_obj formats for data

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_difference_data.m,v 1.4 2007-08-29 09:00:55 aadler Exp $

data_width= max(num_frames(data1), num_frames(data2));

if ~isnumeric(data2) if isfield(data2,'fwd_model')
   if nargin>=3; warning('ignoring parameter 3'); end
   fwd_model= data2.fwd_model;
end ; end
if ~isnumeric(data1) if isfield(data1,'fwd_model')
   if nargin>=3; warning('ignoring parameter 3'); end
   fwd_model= data1.fwd_model;
end ; end

if ~isfield(fwd_model,'normalize_measurements')
   fwd_model.normalize_measurements = 0;
end

fdata1 = filt_data( fwd_model, data1, data_width );
fdata2 = filt_data( fwd_model, data2, data_width );

if fwd_model.normalize_measurements
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end

function nf= num_frames(d0)
   if isnumeric( d0 )
      nf= size(d0,2);
   elseif strcmp( d0(1).type, 'data' );
      nf= length(d0);
   else
      error('Problem calculating number of frames');
   end
   
% test for existance of meas_select and filter data
function d2= filt_data(fwd_model, d0, data_width )
   if ~isnumeric( d0 )
       % we probably have a 'data' object

       l_obj = length(d0);
       d1 = zeros( length( d0(1).meas ), l_obj);
       for i=1:l_obj
          if strcmp( d0(i).type, 'data' )
              d1(:,i) = d0(i).meas;
          else
              error('expecting an object of type data');
          end
       end

   else
      % we have a matrix of data. Hope for the best
      d1 = d0;
   end

   d1= double(d1); % ensure we can do math on our object

   if isfield(fwd_model,'meas_select');
      % we have a meas_select parameter

      meas_select= fwd_model.meas_select;
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

