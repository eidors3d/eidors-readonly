function [elem_data,ref_lev] = scale_for_display( elem_data, ref_lev, clim )
% [elem_data,ref_lev] = scale_for_display( elem_data, ref_lev, clim )
%
% PARAMETERS: elem_data
%  elem_data: data for fem elements or image pixels
%
% PARAMETERS: ref_lev
%  ref_lev:   if param provided, use it,
%               otherwise use the global value
%             Can be numeric or 'auto' or 'use_global' 
%
% PARAMETERS: clim
%    clim - colour limit. Colours more different from ref_level are cropped.
%         - if not specified or scale==[] => no limit
%
% OUTPUT: ref_lev
%    - the actual numerical ref_level used
%
% $Id: scale_for_display.m,v 1.3 2006-08-25 00:11:17 aadler Exp $
% (C) 2006 Andy Adler. Licensed under GPL v2

   global eidors_colours;
   % need crazy if statement to support matlab >= 6.0
   if nargin <=1
      ref_lev = eidors_colours.ref_level;
   elseif strcmp(ref_lev, 'use_global' );
      ref_lev = eidors_colours.ref_level;
   end

   if nargin<=2
      clim= [];
   end

   if ~isnumeric(ref_lev)
      if ~strcmp(ref_lev, 'auto')
          error('ref_level parameter must be "auto" or numeric');
      end
      s_ed= elem_data(:);
      s_ed(isnan(s_ed)) = [];
      s_ed= sort(s_ed);
      e= length(s_ed);
      ref_lev = mean(s_ed( ceil(.25*e):floor(.75*e) ));
   end

   elem_data = elem_data - ref_lev;

   % Crop output to the colour limit
   if ~isempty(clim)
      elem_data( elem_data> clim)=  clim;
      elem_data( elem_data<-clim)= -clim;
   end

