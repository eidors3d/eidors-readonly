function mdl = eidors_model_params( mdl );
% mdl = eidors_model_params( mdl );
% Fill in default parameter values in EIDORS types
%
% (C) 2006 Andy Adler. Licensed under the GPL v2
% $Id: eidors_model_params.m,v 1.1 2006-11-15 19:36:06 aadler Exp $

switch mdl.type;
   case 'fwd_model';
   case 'fwd_model';
      if ~isfield(obj,'normalize_measurements');
         mdl.normalize_measurements= 0;
      end
      mdl.elems=    double(mdl.elems);
      mdl.boundary= double(mdl.boundary);
      % fill in boundary if it doesn't exist
      % fill in boundary n_elems, n_nodes
end
