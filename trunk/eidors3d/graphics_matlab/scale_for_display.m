function elem_data= scale_for_display( elem_data, ref_lev )
% elem_data= scale_for_display( elem_data, ref_lev )
%
%  elem_data: prescaled elem_data
%  ref_lev:   if param provided, use it,
%               otherwise use the global value
%             Can be numeric or 'auto' or 'use_global' 
%
% $Id: scale_for_display.m,v 1.2 2006-08-01 16:33:57 aadler Exp $
% (C) 2006 Andy Adler. Licensed under GPL v2

   global eidors_colours;
   % need crazy if statement to support matlab >= 6.0
   if nargin <=1
      ref_lev = eidors_colours.ref_level;
   elseif strcmp(ref_lev, 'use_global' );
      ref_lev = eidors_colours.ref_level;
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


