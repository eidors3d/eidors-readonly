function   mdl = mash_nodes(mdl, method, idm, dim, elec_true, elec_err)
% function mdl = mash_nodes(mdl, method, idm, dim, elec_true, [elec_err])
%
% Mashes the nodes in the X, Y, or Z direction to
% minimize the error in the electrode locations while
% maintaining the boundaries of a rectangular domain.
%
% This function exists because NetGen and GMSH *suck*
% at getting the electrodes in precisely the correct
% place. In addition it is very expensive to move
% electrodes by a small amount if re-meshing is
% required. Mesh quality should not be adversely
% affected as long as electrodes do not move by too
% much and are not too close to the fixed boundaries.
%
% mdl - a forward or reconstruction (coarse) model
% method - 'shift_all' use when upd_dim == interp_dim
%        to shuffle all nodes along the surface in the
%        direction of the interpolation
%        - 'shift_middle' maintain the outer
%        boundaries and mash the middle around to
%        minimize the error
%        - 'shift_surface' leave the negative boundary
%        alone and shift the positive boundary around
%        to minimize the error
% idm - the dimension in which to index the error (monotonic with electrodes)
% dim - the dimension in which to minimize the error
% elec_true - the true electrode locations
% elec_err - the error in electrode location, single
%        vector (default: error in the upd_dim
%        dimension)
%
% EXAMPLE
% elec_err=1
% while(elec_err > eps)
%  elec_nodes = [fmdl.electrode(:).nodes];
%  elec_err = elec_true - fmdl.nodes(elec_nodes,:);
%  cmdl = mash_nodes(cmdl, 'shift_all',     1, 2, elec_true(:,[2,3]), elec_err(:,[2,3])); % Y (downslope)
%  cmdl = mash_nodes(cmdl, 'shift_surface', 1, 1, elec_true(:,[2,3]), elec_err(:,[2,3])); % Z (vertical)
%  % NOTE that X and Y for cmdl are flipped when plotted in figures but must
%  % remain messed up for dimension order since the coarse2fine mapping function
%  % gives a bad result otherwise
%  fmdl = mash_nodes(fmdl, 'shift_middle',  2, 1, elec_true); % X (cross-slope)
%  fmdl = mash_nodes(fmdl, 'shift_all',     2, 2, elec_true); % Y (downslope)
%  fmdl = mash_nodes(fmdl, 'shift_surface', 2, 3, elec_true); % Z (vertical)
% end
%
% (C) 2014 A. Boyle
  if nargin < 6
     elec_err = mdl_electrode_error(mdl, elec_true);
  end
  err = elec_err(:,dim);

  % add borders for electrode positions at
  % the volume boundary and 50% of the edge to electrode distance
  xq = mdl.nodes(:,idm);
  x = [min(xq); ...
       mean([min(xq) min(elec_true(:,idm))]); ...
       elec_true(:,idm);
       mean([max(xq) max(elec_true(:,idm))]); ...
       max(xq)];
  v = [0; 0; err; 0; 0];
plot(err);
  % scale error to match the electrode locations
  vq = interp1(x, v, xq, 'linear', 'extrap');

  switch method
    case 'shift_all'
      vqs = 1;
    case 'shift_middle'
      yq = mdl.nodes(:,dim);
      yqr = max(yq)-min(yq); % range
      yqm = (max(yq) + min(yq))/2; % middle = (max + min)/2
      vqs = 1-abs(yq-yqm)./(yqr/2); % scale the shift depending on x's distance from midline
    case 'shift_surface' % assumes positive surface
      yq = mdl.nodes(:,dim);
      yqr = max(yq) - min(yq);
      yqm = min(yq); % min
      vqs = abs(yq - yqm)./yqr; % scale by distance from surface
    otherwise
      error(['unrecognized method: ',method]);
  end
  mdl.nodes(:,dim) = mdl.nodes(:,dim) + (vq .* vqs);
end

% calculate the error in electrode position for the fwd_model
function err = mdl_electrode_error(mdl, xyz)
   if ~isfield(mdl, 'electrode')
      error('electrodes not available on this model, must supply positional errror');
   end

   nel=length(mdl.electrode); % number of electrodes
   nd=size(mdl.nodes,2); % number of dimensions

   eu = ones(nel,nd)*NaN; % init
   for i=1:length(mdl.electrode)
      n = mdl.electrode(i).nodes; % nodes per electrode
      eu(i,:) = mean(mdl.nodes(n,:)); % approx centre of each electrode
      % Note: 'eu' could be a *bit* more accurate by dropping the interior nodes on the CEM, but this code is simpler
   end
   err = xyz - eu;
end
