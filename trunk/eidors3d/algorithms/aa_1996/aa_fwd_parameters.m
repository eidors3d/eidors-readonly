function param = aa_fwd_parameters( fwd_model )
% NP_FWD_PARAMETERS: data= np_fwd_solve( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Andy Adler's EIT code
%   param.n_elem   => number of elements
%   param.n_elec   => number of electrodes
%   param.n_node   => number of nodes (vertices)
%   param.n_stim   => number of current stimulation patterns
%   param.n_dims   => dimentions (2= 2D, 3=3D)
%   param.n_meas   => number of measurements (total)
%   param.NODE     => vertex matrix
%   param.ELEM     => connection matrix
%   param.QQ       => Current into each NODE
% $Id: aa_fwd_parameters.m,v 1.3 2005-06-04 17:12:39 aadler Exp $

param = eidors_obj('get-cache', fwd_model, 'aa_1996_fwd_param');

if ~isempty(param)
   eidors_msg('aa_fwd_parameters: using cached value', 3);
   return
end

param = calc_param( fwd_model );

eidors_obj('set-cache', fwd_model, 'aa_1996_fwd_param', param);
eidors_msg('aa_fwd_parameters: setting cached value', 3);

% perform actual parameter calculation
function pp= calc_param( fwd_model );

pp.NODE= fwd_model.nodes';
pp.ELEM= fwd_model.elems';

n= size(pp.NODE,2);        %NODEs
d= size(pp.ELEM,1);        %dimentions+1
e= size(pp.ELEM,2);        %ELEMents
p = length(fwd_model.stimulation );


pp.QQ= zeros(e,p);
for i=1:p
    stim= fwd_model.stimulation(i);
    % for each electrode in this stimulation pattern
    elecs= find( stim.stim_pattern );
    for elec= elecs(:)'
    %   fwd_model.electrode(i).z_contact - unused here
        elec_nodes = fwd_model.electrode(elec).nodes;
        l_elec_nodes = length(elec_nodes);
        curr= stim.stim_pattern( elec );

        pp.QQ( elec_nodes, i) = curr / l_elec_nodes;
    end
end


% pack into a parameter return list
pp.n_elem   = e;
pp.n_elec   = length( fwd_model.electrode );
pp.n_node   = n;
pp.n_stim   = p;
pp.n_dims   = d-1;
%param.n_meas   = 
