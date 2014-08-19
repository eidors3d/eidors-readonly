function   mdl = mash_nodes(mdl, method, interp_dim, upd_dim, elec_true, elec_err)
% function mdl = mash_nodes(mdl, method, interp_dim, upd_dim, elec_true, [elec_err])
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
% interp_dim - the dimension in which to interpolate
% upd_dim - the dimension in which to minimize the
%        error
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
%  cmdl = mash_nodes(cmdl, 'shift_all',     2, elec_true, elec_err(:,2)); % Y (downslope)
%  cmdl = mash_nodes(cmdl, 'shift_surface', 1, elec_true, elec_err(:,3)); % Z (vertical)
%  % NOTE that X and Y for cmdl are flipped when plotted in figures but must
%  % remain messed up for dimension order since the coarse2fine mapping function
%  % gives a bad result otherwise
%  fmdl = mash_nodes(fmdl, 'shift_middle',  1, elec_true, elec_err(:,1)); % X (cross-slope)
%  fmdl = mash_nodes(fmdl, 'shift_all',     2, elec_true, elec_err(:,2)); % Y (downslope)
%  fmdl = mash_nodes(fmdl, 'shift_surface', 3, elec_true, elec_err(:,3)); % Z (vertical)
% end
%
% (C) 2014 A. Boyle
  if nargin < 5
     if isfield(mdl, 'electrode')
        error('electrodes not available on this model, must supply positional errror');
     else
        % this code assumes a single node per electrode --- Point Electrode Model
        % MATLAB voodoo: convert struct array to array
        elec_nodes = [mdl.electrode(:).nodes];
        elec_err = elec_true - mdl.nodes(elec_nodes,:);
        elec_err = elec_err(:,upd_dim);
     end
  end

  y = mdl.nodes(:,interp_dim);
  yi = [min(y); ...
        mean([min(y) min(elec_true(:,2))]); ...
        elec_true(:,interp_dim);
        mean([max(y) max(elec_true(:,2))]); ...
        max(y)];
  yv = [0; 0; elec_err; 0; 0];

  % scale error to match the electrode locations
  yy = interp1(yi, yv, y, 'linear', 'extrap');

  switch method
    case 'shift_all'
      mdl.nodes(:,upd_dim) = mdl.nodes(:,upd_dim) + yy;
    case 'shift_middle'
      x = mdl.nodes(:,upd_dim);
      xr = range(x);
      xm = xr/2 + min(x); % middle
      xs = 1-(x-xm)./(xr/2); % scale the shift depending on x's distance from midline
      mdl.nodes(:,upd_dim) = x + yy.*xs;
    case 'shift_surface' % assumes positive surface
      z = mdl.nodes(:,upd_dim);
      zr = range(z);
      zm = min(z); % min
      zs = (z - zm)/zr; % scale by distance from surface
      mdl.nodes(:,upd_dim) = z + yy.*zs;
    otherwise
      error(['unrecognized method: ',method]);
  end
end

