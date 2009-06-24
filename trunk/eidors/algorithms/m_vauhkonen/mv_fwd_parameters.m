function p= mv_fwd_parameters( fmdl )
% MV_FWD_PARAMETERS: data= np_fwd_solve( fwd_model )
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Marco Vauhkonen Polydorides EIDORS2D code
%   param.n_elem   => number of elements
%   param.n_elec   => number of electrodes
%   param.n_node   => number of nodes (vertices)
%   param.n_stim   => number of current stimulation patterns
%   param.n_meas   => number of measurements (total)
%   param.vtx      => vertex matrix
%   param.simp     => connection matrix
%   param.srf      => boundary triangles
%   param.df       => vector of measurements for each current pattern
%   param.elec     => nodes attached to each electrode
%   param.zc       => vector of contact impedances
%   param.indH     => electrodes used for each measurement
%   param.I        => RHS (current term) for FEM solution
%   param.Ib       => Current for electrodes
%   param.perm_sym => 'sym' parameter
%   param.gnd_ind  => node attached to ground
%   param.normalize  => difference measurements normalized?

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

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
       


function eno=in_electrode(face, fmdl)
   eno= 0;
   for i=1:length(fmdl.electrode)
      is = intersect(face(:), fmdl.electrode(i).nodes(:) );
      if length(is)==2; eno=i; end
   end
     
