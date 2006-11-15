function mdl = eidors_model_params( mdl );
% mdl = eidors_model_params( mdl );
% Fill in default parameter values in EIDORS types
%
% (C) 2006 Andy Adler. Licensed under the GPL v2
% $Id: eidors_model_params.m,v 1.3 2006-11-15 19:45:33 aadler Exp $

% TODO - to caching

switch mdl.type;
   case 'inv_model';
      mdl.fwd_model= eidors_model_params( mdl.fwd_model );

   case 'fwd_model';
      if ~isfield(mdl,'normalize_measurements');
         mdl.normalize_measurements= 0;
      end
      mdl.elems=    double(mdl.elems);
      mdl.boundary= double(mdl.boundary);

      % fill in boundary if it doesn't exist
      mdl.n_elem = size(mdl.elems,1);
      mdl.n_node = size(mdl.nodes,1);
      mdl.n_elec = length(mdl.electrode);
end
