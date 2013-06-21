function data = fwd_solve_elec_move(jnk, img)

img = convert_legacy_movement_img(img);

if isfield(img, 'movement')
   img = apply_movement(img);
   img = rmfield(img,'movement');
end
try
   solver = img.fwd_model.fwd_solve_elec_move.data_solver;
catch
   solver = 'eidors_default';
end
img.fwd_model.solve = solver;
data = fwd_solve(img);
data.name = 'Solved by fwd_solve_elec_move';
function img = apply_movement(img)

n_elec  = length(img.fwd_model.electrode);
n_nodes = length(img.fwd_model.nodes);
n_dim   = mdl_dim(img.fwd_model);

try 
   node_move = img.movement.node_data;
catch
   node_move = zeros(n_nodes, n_dim );
end

if isfield(img.movement, 'electrode_data');
   move = reshape( img.movement.electrode_data , n_elec, n_dim );
   for i = 1:n_elec
      enodes = img.fwd_model.electrode(i).nodes;
      lenodes = length(enodes);
      node_move(enodes,:) = node_move(enodes,:) + repmat(move(i,:),lenodes,1);
   end
end

img.fwd_model.nodes = img.fwd_model.nodes + node_move;
