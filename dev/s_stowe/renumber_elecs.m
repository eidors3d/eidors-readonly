function [fmdl] = renumber_electrodes( varargin ) 
% Renumbers the electrodes of a forward model based on a given list of electrode centres or
% Attempts to automatically number electrodes based on selected pattern.
% Usage
% fwd_model = renumber_electrodes(fwd_model, electrode_coordinates, coordinate_format)
% Give the coordinates in cylindrical or Cartesian coordinates and select format
% electrode_coordinates -> (r, angle, z; r, angle, z; ...) or (x, y, z; x, y, z; ...)
% coordinate_format -> (0 for Cartesian, 1 for Cylindrical)
% fwd_model = renumber_electrodes(fwd_model, opts)
% Attempts to renumber electrodes automatically based on the following user options:
% Renumber options
% opts.ext_elecs -> number of external electrodes
% opts.ext_rows -> number of vertical electrode rows
% opts.int_elecs -> number of internal electrodes
% opts.offset -> electrode offset from centre (sternum) in degrees (default 0)
% opts.pattern -> application pattern for external electrodes (default 'square')
	if (length(varargin) == 3)  % Use specified coordinates
		fmdl = varargin{1};
		coordinates = varargin{2};
		coordinate_type = varargin{3};
		old_electrode = fmdl.electrode; % Keep original order
		if isempty(horzcat(old_electrode.nodes))
		    for i=1:size(old_electrode,2)
		        old_electrode(i).nodes = unique(old_electrode(i).faces,'stable')';
		    end
		end
		elec_nodes = horzcat(old_electrode.nodes);
		elec_ind = [];
		for i=1:size(old_electrode,2)
		    elec_ind = [elec_ind; ones(size(old_electrode(i).nodes,2),1)*i];
		end
		
		fmdl.electrode = [];
		% Convert all nodes to coordinates
		elec_coords = zeros(length(elec_nodes),3);
		for i=1:length(elec_nodes)
		    elec_coords(i,:) = fmdl.nodes(elec_nodes(i),:);
		end
		elec_coords = [elec_ind elec_coords];
		used_elecs = [];
	    for i=1:size(coordinates,1)
			% Find the nearest electrode node and assign the label of current coordinates
			if(coordinate_type == 0) % Cartesian
				approx_node(1,:) = [0 coordinates(i,:)];
			elseif(coordinate_type == 1) % Cylindrical
				approx_node(1,1) = 0;
				approx_node(1,2) = coordinates(i,3)*cosd(coordinates(i,2)); % x location
				approx_node(1,3) = coordinates(i,3)*sind(coordinates(i,2)); % y location
				approx_node(1,4) = coordinates(i,3);
			else
				error('Only coordinate types 0 (Cartesian) or 1 (cylindrical) are supported.');
			end
	        compare = abs(elec_coords-approx_node);
	        sum_compare = sum(compare(:,2:4),2);
	        [loc,~] = find(sum_compare == min(sum_compare));
	        closest_elec = elec_ind(loc);
	        if isempty(horzcat(old_electrode.nodes))
	            fmdl.electrode(i).faces = old_electrode(closest_elec).faces;
	            fmdl.electrode(i).nodes = [];
	        else
	            fmdl.electrode(i).nodes = old_electrode(closest_elec).nodes;
	        end
	        fmdl.electrode(i).z_contact = old_electrode(closest_elec).z_contact;
	        used_elecs = [used_elecs unique(closest_elec)];
	    end
	elseif (length(varargin) == 2)  % Use specified options
		fmdl = varargin(1);
		opts = varargin(2);
		old_electrode = fmdl.electrode; % Keep original order
		if isempty(horzcat(old_electrode.nodes))
		    for i=1:size(old_electrode,2)
		        old_electrode(i).nodes = unique(old_electrode(i).faces,'stable')';
		    end
		end
		elec_nodes = horzcat(old_electrode.nodes);
		elec_ind = [];
		for i=1:size(old_electrode,2)
		    elec_ind = [elec_ind; ones(size(old_electrode(i).nodes,2),1)*i];
		end
		
		fmdl.electrode = [];
		% Convert all nodes to coordinates
		elec_coords = zeros(length(elec_nodes),3);
		for i=1:length(elec_nodes)
		    elec_coords(i,:) = fmdl.nodes(elec_nodes(i),:);
		end
		elec_coords = [elec_ind elec_coords];
		
		% Get all the angles
		divs = opts.ext_elecs/opts.ext_rows;
		sep_angle = 360/divs;  % TODO - Insert check for even division
		node_compare = zeros(size(elec_coords));
		used_elecs = [];
		if opts.pattern == 'square'
		    for i=1:ext_elecs
		        angle = (ceil(i/2)-1)*(sep_angle+opts.offset);
		        ud_loc = (ceil((i+1)/2-1));
		        if ~mod(ud_loc,2)
		            height = max(elec_coords(:,4));
		        end
		        if mod(ud_loc,2)
		            height = min(elec_coords(:,4));
		        end
		        node_compare(:,2) =  max(elec_coords(:,2))*sind(angle); % x location
		        node_compare(:,3) = -max(elec_coords(:,3))*cosd(angle); % y location
		        node_compare(:,4) = height; % z location
		        compare = abs(elec_coords-node_compare);
		        sum_compare = sum(compare(:,2:4),2);
		        [loc,~] = find(sum_compare == min(sum_compare));
		        closest_elec = elec_ind(loc);
		        if isempty(horzcat(old_electrode.nodes))
		            fmdl.electrode(i).faces = old_electrode(closest_elec).faces;
		            fmdl.electrode(i).nodes = [];
		        else
		            fmdl.electrode(i).nodes = old_electrode(closest_elec).nodes;
		        end
		        fmdl.electrode(i).z_contact = old_electrode(closest_elec).z_contact;
		        used_elecs = [used_elecs closest_elec];
		    end
		else
			error('This electrode pattern is not yet supported.');
		end
		
		subplot(1,2,2)
		show_fem(fmdl,[0 1.02])
		
		used_idx = find(~ismember(elec_coords(:,1),used_elecs));
		remaining_elec_nodes = elec_coords(used_idx,:);
		remaining_elecs = unique(remaining_elec_nodes(:,1));
		
		for i=ext_elecs+1:int_elecs+ext_elecs
		    [loc,~] = find(remaining_elec_nodes(:,4) == max(remaining_elec_nodes(:,4)));
		    highest_elec = remaining_elec_nodes(loc,1);
		    if length(unique(highest_elec)) ~= 1;
				keyboard
		    end
		    fmdl.electrode(i).nodes = old_electrode(highest_elec).nodes;
		    fmdl.electrode(i).z_contact = old_electrode(highest_elec).z_contact;
		    [temp,~] = find(remaining_elec_nodes(:,1) == unique(highest_elec));
		    remaining_elec_nodes(temp,:) = [];
		end

	else
		error('Renumber electrodes requires 2 or 3 inputs. See "help renumber_electrodes" for details.')
	end
end

