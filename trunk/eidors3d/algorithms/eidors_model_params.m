function mdl = eidors_model_params( mdl );
% mdl = eidors_model_params( mdl );
% Fill in default parameter values in EIDORS types
%
% (C) 2006 Andy Adler. Licensed under the GPL v2
% $Id: eidors_model_params.m,v 1.8 2008-07-10 15:10:39 aadler Exp $

% TODO - to caching

try
   type=mdl.type;
catch
   error('eidors_model_params: object is not eidors object (no type)');
end

switch type;
   case 'inv_model';
      mdl.fwd_model= eidors_model_params( mdl.fwd_model );
      if ~isfield(mdl,'reconst_type');
         mdl.reconst_type= 'difference';
      end

   case 'fwd_model';
      if ~isfield(mdl,'normalize_measurements');
         mdl.normalize_measurements= 0;
      end
      if ~isfield(mdl,'elems');
         return;
      end
      
      mdl.elems=    double(mdl.elems);
      if isfield(mdl,'boundary')
         mdl.boundary= double(mdl.boundary);
      else
         mdl.boundary= find_boundary(mdl.elems);
      end

      % fill in boundary if it doesn't exist
      mdl.n_elem = size(mdl.elems,1);
      mdl.n_node = size(mdl.nodes,1);
      if isfield(mdl,'electrode');
         mdl.n_elec = length(mdl.electrode);
      else
         mdl.n_elec = 0;
      end
end
