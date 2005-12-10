function img = inv_solve( inv_model, data1, data2)
% INV_SOLVE: calculate imag from an inv_model and data
% 
% inv_solve can be called as
%     img= inv_solve( inv_model, data1, data2)
%   if inv_model.type = 'difference'
% or
%     img= inv_solve( inv_model, data )
%   if inv_model.type = 'static'
%
% in each case it will call the inv_model.solve
%
% data      is a measurement data structure
% inv_model is a inv_model structure
% img       is an image structure
%           or a vector of images is data1 or data2 are vectors
%
% For difference EIT:
% data1      => difference data at earlier time
% data2      => difference data at later time
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

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: inv_solve.m,v 1.13 2005-12-10 10:47:29 aadler Exp $

% COMMENT: There seems to be no general way to cache
%       inv_model parameters. Thus, each algorithm needs
%       to do it separately. It may be possible to cache
%       one-step inverse matrices, but that is not done. 

eidors_msg(['inv_solve:', inv_model.name],1);

% caching of images is disabled -
%  in order to implement image caching, it needs to be possible
%  to tag them with both data1 and data2


if     strcmp(inv_model.reconst_type,'static')
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
   fdata1 = filt_data( inv_model, data1 ); l_data1= size(fdata1,2);
   fdata2 = filt_data( inv_model, data2 ); l_data2= size(fdata2,2);

   if l_data1 ~=1 & l_data2 ~=1 & l_data1 ~= l_data2
      error('inconsistent number of specified measurements');
   elseif l_data1 >1 & l_data2==1 
      fdata2 = fdata2 * ones(1,l_data1);
   elseif l_data1==1 & l_data2 >1 
      fdata1 = fdata1 * ones(1,l_data2);
   end

   imgc= feval( inv_model.solve, inv_model, fdata1, fdata2);
else
   error('inv_model.reconst_type not understood'); 
end

elem_data= imgc.elem_data;
n_img = size(elem_data,2);

for i=1:n_img
   imgc.elem_data= elem_data(:,i);
   img(i) = eidors_obj('image', imgc );
end

% test for existance of meas_select and filter data
function d1= filt_data(inv_model, d0 )
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

    elseif isfield(inv_model.fwd_model,'meas_select');
       % we have a meas_select

       meas_select= inv_model.fwd_model.meas_select;
       if     size(d0,1) == length(meas_select)
          d1= d0(meas_select,:);
       elseif size(d0,1) == sum(meas_select>0)
          d1= d0;
       else
          error('data size does not match meas_select');
       end

    else
       % we have no info about whether data is ok
       % hope for the best
       d1 = d0;
    end
    d1= double(d1); % ensure we can do math on our object
