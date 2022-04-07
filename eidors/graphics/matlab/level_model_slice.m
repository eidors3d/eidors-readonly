function [out, out2, out3] = level_model_slice(varargin)
%LEVEL_MODEL_SLICE - level 3D points for slicing at z=0
% 
% XYZ = LEVEL_MODEL_SLICE(NODES, LEVEL)
% XYZ = LEVEL_MODEL_SLICE(NODES, LEVEL, SN)
% Transforms a list of 3D points such that level is a slice at z=0.
%	NODES: n_nodes x 3 list of point coordinates 
%	LEVEL: definition of the slice(s). Can be:
%     - scalar: number of slices along the z axis
%     - vector (N x 3): intercepts with the coordinate axis
%     - struct:
%       .centre (N x 3): centre point for transformation
%		.normal_angle (1 x 4) [nx, ny, nz, theta]: normal vector and angle 
%		  of rotation around the axis defined by it (right hand rule).
%         The model will be rotated such that the normal vector becomes 
%         [0,0,1].         
%		.rotation_matrix (3 x 3)
%		An error is thrown if both normal_angle and rotation_matrix are present.
%   SN	 : slice number to use if level defines multiple slices (default: 1).
%	XYZ	 : NODES transformed such that level is a slice at z=0.
%
% MDL = LEVEL_MODEL_SLICE(MDL,__) 
% Modifies the nodes of an EIDORS fwd_model structure MDL. Can be used in
% place of NODES with every other syntax.
%
% N = LEVEL_MODEL_SLICE(LEVEL) returns the number of slices defined in LEVEL.
%
% [C,M] = LEVEL_MODEL_SLICE(LEVEL, SN) (discouraged) returns the center C and
% rotation matrix M of slice SN (optional, default 1). An error is thrown if
% level is scalar. Note that two outputs must always be requested for this 
% syntax.  
%
% [C,M] = LEVEL_MODEL_SLICE(NODES, LEVEL, SN) (discouraged) allows to obtain C
% and M when LEVEL is scalar. Note that two outputs must always be requested for
% this syntax. 

% (C) 2006-2021 Andy Adler and Bartek Grychtol 
% License: GPL version 2 or version 3
% $Id$

if nargin==1 && isstr(varargin{1}) && strcmp(varargin{1}, 'UNIT_TEST')
  do_unit_test; return
end

if nargout == 1 || nargout == 3
	if nargin == 1
        % N = LEVEL_MODEL_SLICE(LEVEL) returns the number of slices defined in LEVEL
        out = get_number_of_slices(varargin{1});
	else
		% XYZ = LEVEL_MODEL_SLICE(NODES, LEVEL)
		% XYZ = LEVEL_MODEL_SLICE(NODES, LEVEL, SN)
        [C, M] = matrix_from_level(varargin{:});
        SHIFT = [C(1), C(2), C(3)] * M';
        SHIFT(3) = 0;
        if isstruct(varargin{1})
            % accept fwd_model as input
            out = varargin{1};
            % (nodes - C) * M' + SHIFT
            out.nodes = bsxfun(@plus,bsxfun(@minus,out.nodes, C ) * M', SHIFT);
        else
            % (nodes - C) * M' + SHIFT
            out = bsxfun(@plus, bsxfun(@minus, varargin{1}, C ) * M', SHIFT);
        end
        if nargout == 3
            out2 = C;
            out3 = M;
        end
        
	end
elseif nargout == 2
	if nargin < 3
		% [C,M] = LEVEL_MODEL_SLICE(LEVEL)
		% [C,M] = LEVEL_MODEL_SLICE(LEVEL, SN)
        [out, out2] = matrix_from_level(varargin{:});
	else
		% [C,M] = LEVEL_MODEL_SLICE(NODES, LEVEL, SN)
        [out, out2] = matrix_from_number(varargin{:});
	end
end

function N = get_number_of_slices(level)  
	if isnumeric(level) && isscalar(level)
		N = level;
	elseif isnumeric(level)
		N = size(level,1);
	elseif isstruct(level) && isscalar(level)
		N = size(level.centre, 1);
	else
		error('level must be numeric or a scalar struct');
	end

function [C,M] = matrix_from_level(varargin)
    nodes = []; sn = 1;
    switch nargin
        case 3 % nodes, level, sn
            nodes = varargin{1}; level = varargin{2}; sn = varargin{3};
        case 2 % level & sn or nodes & level
            if isnumeric(varargin{2}) && isscalar(varargin{2}) 
                level = varargin{1}; sn = varargin{2};
            else
                nodes = varargin{1}; level = varargin{2};
            end
        case 1 % level
            level = varargin{1};
    end

    if isnumeric(level) && isscalar(level)
        if isempty(nodes)
            error('Need nodes to calculate with scalar numeric level');
        else
            [C, M] = matrix_from_number(nodes, level, sn);
        end
    elseif isnumeric(level)
        [C, M] = matrix_from_intercept(level, sn);
    elseif isstruct(level) && isscalar(level)
        [C, M] = matrix_from_struct(level, sn);
    else
        error('level must be numeric or a scalar struct');
    end

    
function [C, M] = matrix_from_number(nodes, num_levels, sn)
    if isstruct(nodes)
        nodes = nodes.nodes;
    end
    M = diag([1,1,1]);
    
    zmax= max(nodes(:,3));
    zmin= min(nodes(:,3));
    levels = linspace(zmin,zmax, num_levels+2);
    levels = levels(2:end-1);
    
    C = (max(nodes) + min(nodes))/2;
    C(3) = levels(sn);
    
    
    
  
function [C, M] = matrix_from_intercept(level, sn)   
   if nargin < 2
     sn = 1;
   end
   level = level(sn, :);
   R = eye(3);
   if length(level)==4
     theta = -level(4); % we rotate the model, user thinks of the image
     R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0,0,1;];
   end
   
   % Infinities tend to cause issues -> replace with realmax
   % Don't need to worry about the sign of the inf
   level( isinf(level) | isnan(level) ) = realmax;
   level( level==0 ) =     1e-10; %eps;

   % Step 1: Choose a centre point in the plane
   %  Weight the point by it's inv axis coords
   invlev= 1./level;
   ctr= invlev / sum( invlev.^2 );

   % Step 2: Choose basis vectors in the plane
   %  First is the axis furthest from ctr
   [jnk, s_ax]= sort( - abs(level - ctr) );
   v1= [0,0,0]; v1(s_ax(1))= level(s_ax(1));
   v1= v1 - ctr;
   v1= v1 / norm(v1);

   % Step 3: Get off-plane vector, by cross product
   v2= [0,0,0]; v2(s_ax(2))= level(s_ax(2));
   v2= v2 - ctr;
   v2= v2 / norm(v2);
   v3= cross(v1,v2);

   % Step 4: Get orthonormal basis. Replace v2
   v2= cross(v1,v3);

   % Step 5: Get bases to point in 'positive directions'
   v1= v1 * (1-2*(sum(v1)<0));
   v2= v2 * (1-2*(sum(v2)<0));
   v3= v3 * (1-2*(sum(v3)<0));
   
   C = ctr;
   M = [v1;v2;v3];
%   NODE= [v1;v2;v3] * (vtx' - ctr'*ones(1,nn) );
  
    
function [C, M] = matrix_from_struct(level, sn)
  if nargin < 2
    sn = 1;
  end
  if ~isfield(level,'centre')
      error('Struct input ''level'' not understood.');
  end
  C = level.centre(sn,:); C = C(:)'; % make row
  
  if isfield(level, 'rotation_matrix')
      M = level.rotation_matrix;
      return % no questions asked
  end
  
  if isfield(level, 'normal') && ~isfield(level,'normal_angle')
      level.normal_angle = level.normal; %quietly accept 
  end
 
  if ~isfield(level,'normal_angle')
      error('Struct input ''level'' not understood.');
  end
  
  A = level.normal_angle(1:3);
  normA = norm(A);
  if normA == 0
    error('Normal vector has zero norm')
  end
  A = A(:) / normA; % normalize and make column
  B = [0,0,1]';
  R = eye(3);
  if length(level.normal_angle) == 4
    theta = level.normal_angle(4);
    R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0,0,1;];
  end
   
  cross_AB = cross(A,B);
  if norm(cross_AB) == 0 % parallel vectors
    if A(3) > 0
      M = R*eye(3);
    elseif A(3) < 0
      M = R*diag([1,-1,-1]);
    end
    return
  end
  
  if 0
      % from https://math.stackexchange.com/questions/180418/
      % calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
      % works well, but leaves a rotation that's hard to control
      dot_AB = A(3);
      ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
      M = eye(3) + ssc(cross_AB) + ...
          ssc(cross_AB)^2*(1-dot_AB)/(norm(cross_AB)^2);
  else
      % here we define a nice coordinate system on the plane, and then
      % transform it to become the identity (ie. invert).
      
      % project each coordinate axis on the plane
      I = eye(3);
      for i = 1:3
          proj(:,i) = I(:,i) - (dot(I(:,i),A) / normA^2) * A;
      end
      norm_proj = vecnorm(proj);
      max_norm = max(norm_proj);
      
      % choose what happens based on which projection is longest
      M = zeros(3);
      M(:,3) = A;
      if norm_proj(3) == max_norm
          % projection of z becomes y
          M(:,2) = proj(:,3);
          M(:,1) = cross(M(:,2),M(:,3));
      elseif norm_proj(2) == max_norm
          % projection of y becomes y
          M(:,2) = proj(:,2);
          M(:,1) = cross(M(:,2),M(:,3));
      else
          % projection of x becomes x
          M(:,1) = proj(:,1);
          M(:,2) = cross(M(:,3),M(:,1));
      end
      M = inv(M);  
  end
  M = R*M ;  
  

  
function do_unit_test
  unit_test_cmp('Scalar level', level_model_slice(5),5);
  unit_test_cmp('Vector level', level_model_slice([inf, inf, 3]),1);
  unit_test_cmp('Matrix level', level_model_slice([inf, inf, 3; inf, inf, 5]),2);
  LVL = [inf, inf, 3; inf, inf, 5];
  LVL(:,4) = rand(1,2);
  unit_test_cmp('Matrix level angle', level_model_slice(LVL),2);
  lvl.centre = [ 0,0,0];
  unit_test_cmp('Struct level', level_model_slice(lvl),1);
  lvl.centre(2,:) = [0,0,3];
  unit_test_cmp('Struct level', level_model_slice(lvl),2);
  lvl.normal_angle = rand(3,1)-.5;
  lvl.normal_angle = lvl.normal_angle / norm(lvl.normal_angle);
  lvl.normal_angle(4) = 2*pi * rand(1);
  [C, M] = level_model_slice(lvl);
  unit_test_cmp('Normal -> Matrix', M*lvl.normal_angle(1:3),[0,0,1]', 5*eps);
  [C, M] = level_model_slice([inf 5 inf]);
  lvl.normal_angle = [0,0,-1];
  [C, M] = level_model_slice(lvl);
  unit_test_cmp('Parallel normal', M,diag([1,-1,-1]));
  
  lvl.normal_angle = [0,0,0];
  try
    [C, M] = level_model_slice(lvl);
    fprintf('TEST: %20s = %4s \n','Expect error', 'Fail');  
  catch
    fprintf('TEST: %20s = %4s \n','Expect error', 'OK'); 
  end
