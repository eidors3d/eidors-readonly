function param= mv_fwd_parameters( fwd_model )
% MV_FWD_PARAMETERS: data= mv_fwd_solve( fwd_model )
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Marco Vauhkonen Polydorides EIDORS2D code
%   param.n_elem   => number of elements
%   param.n_elec   => number of electrodes
%   param.n_node   => number of nodes (vertices)
%   param.n_stim   => number of current stimulation patterns
%   param.Element  => Element structure (with Topology and Faces)
%   param.z_contact=> column vector of contact impedances

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$


param = eidors_obj('get-cache', fwd_model, 'mv_fwd_parameters');

if ~isempty(param)
   eidors_msg('mv_fwd_parameters: using cached value', 3);
   return
end

param = calc_param( fwd_model );

eidors_obj('set-cache', fwd_model, 'mv_fwd_parameters', param);
eidors_msg('mv_fwd_parameters: setting cached value', 3);

function p= calc_param( fmdl );

p.n_elem = size(fmdl.elems,1);
p.n_elec = length(fmdl.electrode);
p.n_node = size(fmdl.nodes,1);
p.n_stim = length(fmdl.stimulation);
for i=1:p.n_elem
   elem = fmdl.elems(i,:);
   p.Element(i).Topology = elem;
   p.Element(i).Face =  cell(3,3); % Middle column isn't used
   j=1;for choose = [[1,2];[2,3];[3,1]]'
      ch_face = elem(choose);
      p.Element(i).Face{j,1} = ch_face;
      p.Element(i).Face{j,3} = in_electrode(ch_face, fmdl);
      j=j+1;
   end
end
       
for i=1:p.n_node
   nn.Coordinate = fmdl.nodes(i,:);
   nn.ElementConnection = find(any( fmdl.elems == i, 2))';
   Nodes = fmdl.elems( nn.ElementConnection,: );
   Nodes = unique(Nodes(:))';
   Nodes( Nodes == i ) = [];
   nn.NodeConnection = Nodes;
   p.Node(i)= nn;
end

p.z_contact = vertcat([fmdl.electrode(:).z_contact]);


function eno=in_electrode(face, fmdl)
   eno= 0;
   for i=1:length(fmdl.electrode)
      is = intersect(face(:), fmdl.electrode(i).nodes(:) );
      if length(is)==2; eno=i; end
   end
     
