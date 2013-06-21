function img = convert_legacy_movement_img(img)

% do nothing if the image doesn't need conversion
% old movement images all had elem_data
if ~isfield(img,'elem_data'), return;  end
% check if elem_data is longer than usual
exp_sizes = expected_data_sizes(img);
n_elem = size(img.elem_data,1);
if any( exp_sizes == n_elem), return, end;

n_data_elem = exp_sizes( find(exp_sizes < n_elem, 1,'last') );
elem_data   = img.elem_data;

img = rmfield(img, 'elem_data');
img.conductivity.elem_data  = elem_data(1:n_data_elem,:);
img.movement.electrode_data = elem_data((n_data_elem+1):end, :);
img.physics_data_mapper = 'conductivity';
img.fwd_model.solve     = @fwd_solve_elec_move;

function exp_sizes = expected_data_sizes(img)
% number of elements
exp_sizes(1) = size(img.fwd_model.elems);
% number of elements in coarse model
try
   exp_sizes(2) = size(img.fwd_model.coarse2fine,2);
end

exp_sizes = sort(exp_sizes);