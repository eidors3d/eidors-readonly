% Module: CreateTankWithCompoundElectrodes.m
%
% Creates cylindrical tank phantom with compound electrodes 
%
% Known issues / Todo: -
%
% Input:    Tank - structure
%               .Height
%               .Radius
%           	.maxMesh - 0 for turning off refinement
%
%           Electrodes - structure
%               .InnerRadius1   -   Inner radius of the inner disk
%               .InnerRadius2   -   Outer radius of the inner disk
%               .OuterRadius1   -   Inner radius of the outer disk
%               .OuterRadius2   -   Outer radius of the outer disk
%               .maxh - 0 for turning off refinement
%
% Return:   fmdl
%
% Written by Steffen Kaufmann <sk@steffen-kaufmann.com>, mainly based on 
% http://eidors3d.sourceforge.net/tutorial/netgen/netgen_gen_models.shtml
% Demo #17 from Andy Adler
% March 2014
function fmdl = CreateTankWithCompoundElectrodes(Tank, Electrodes)
    % Calculate Electrode positions
    elec_pos = CalculateElectrodePositions(Electrodes.NumberOf, Electrodes.ZPositions, Tank.Radius);
    
    % Create Netgen string for creating the tank phantom with the compound electrodes
    shape_str = [ng_make_Cylinder(Tank.Height, Tank.Radius, Tank.maxMesh) ng_make_CompoundElectrodes(elec_pos, Electrodes.InnerRadius1, Electrodes.InnerRadius2, Electrodes.OuterRadius1, Electrodes.OuterRadius2, Electrodes.maxh)];

    % A compound electrode consits of two electrodes, therefore we have to
    % duplicate elec_pos - first colum points to the FEM - seem to be
    % needed see demo
    elec_pos = [elec_pos(1, :); elec_pos(2:end, :); elec_pos(2:end, :)];

    % Std. settings according to demo
    elec_shape=[0.1];
    elec_obj = 'cyl';
    
    % Generate model based on previous settings
    fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
    
     % Throw away the first electrode, according to the demo
    fmdl.electrode = fmdl.electrode(2:end);
end

% Calculates the electrode positions in 3d space (first colum points to the 
% FEM - seem to be needed see demo)
function elec_pos = CalculateElectrodePositions(n_elec, Z_position, Radius)
    th = linspace(0,2*pi, n_elec+1)'; th(end)=[];
    on_elecs = ones(n_elec, 1);
    el_th = []; 
    el_z  = []; 
    el_th = [el_th; th];
    el_z  = [el_z ; on_elecs*Z_position];
    
    elec_pos = [Radius*cos(el_th),Radius*sin(el_th),el_z, ones(n_elec, 1)*NaN, ones(n_elec, 1)*NaN, ones(n_elec, 1)*NaN];    
    elec_pos = [0, -Radius, 0, NaN, NaN, NaN; elec_pos];    % get rest of tank first, to remove
end

% Creates the netgen string for compound electrodes
function CMD=ng_make_CompoundElectrodes(elec_pos, InnerRadius1, InnerRadius2, OuterRadius1, OuterRadius2, maxh)
    CMD = '';
    
   if maxh==0; 
      maxh = '';
   else
      maxh = sprintf('-maxh=%f',maxh);
   end    
   
    for electrode=1:size(elec_pos,1)
        x = elec_pos(electrode,1);
        y = elec_pos(electrode,2);
        z = elec_pos(electrode,3);
        
        name = sprintf('elec%04d',electrode);

        if (InnerRadius1 > 0)
            CMD = [CMD sprintf(['solid in_elec = sphere(%6.3f,%6.3f,%6.3f;%6.3f)' ...
                                'and not sphere(%6.3f,%6.3f,%6.3f;%6.3f) %s;\n' ...
                                'solid in_%s= in_elec and mainobj;\n' ...
                                'tlo in_%s cyl;\n' ...
                                'solid out_elec = sphere(%6.3f,%6.3f,%6.3f;%6.3f)' ...
                                'and not    sphere(%6.3f,%6.3f,%6.3f;%6.3f) %s;\n' ...
                                'solid out_%s= out_elec  and mainobj;\n' ...
                                'tlo out_%s cyl;\n'], x,y,z, InnerRadius2, x,y,z, InnerRadius1, maxh, name, name, x,y,z, OuterRadius2, x,y,z, OuterRadius1, maxh, name,name)];
        else
            CMD = [CMD sprintf(['solid in_elec = sphere(%6.3f,%6.3f,%6.3f;%6.3f)' ...
                                ' %s;\n' ...
                                'solid in_%s= in_elec and mainobj;\n' ...
                                'tlo in_%s cyl;\n' ...
                                'solid out_elec = sphere(%6.3f,%6.3f,%6.3f;%6.3f)' ...
                                'and not    sphere(%6.3f,%6.3f,%6.3f;%6.3f) %s;\n' ...
                                'solid out_%s= out_elec  and mainobj;\n' ...
                                'tlo out_%s cyl;\n'],x,y,z, InnerRadius2, maxh, name, name, x,y,z,OuterRadius2, x,y,z, OuterRadius1, maxh, name, name)];
        end;
    end;
end

% Creates netgen string for the tank
function CMD=ng_make_Cylinder(tank_height, tank_radius, maxsz)
   if maxsz==0; 
      maxsz = '';
   else
      maxsz = sprintf('-maxh=%f',maxsz);
   end

   CMD = sprintf(['solid cyl=cylinder (0,0,0;0,0,%6.3f;%6.3f); \n', ...
            'solid mainobj= plane(0,0,0;0,0,-1)\n' ...
                'and  plane(0,0,%6.3f;0,0,1)\n' ...
                'and  cyl %s;\n'],tank_height, tank_radius, tank_height, maxsz);  

end
