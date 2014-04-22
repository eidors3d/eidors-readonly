el_pos = [190,0.5;170,0.5];
extra = {'cube',['solid cube = orthobrick(-0.08,-0.98,0.6 ;0.17,0,0.7) or ' ...
                              'orthobrick( 0.14,-0.98,0.45;0.17,0,0.7) or ' ...
                              'orthobrick(-0.08,-0.98,0.45;-0.05,0,0.7)  or ' ...
                              'orthobrick(-0.17,-0.98,0.3 ;-0.14,0,0.55) or ' ...
                              'orthobrick( 0.05,-0.98,0.3 ;0.08,0,0.55)  or ' ...
                              'orthobrick(-0.17,-0.98,0.3 ;0.08,0,0.4);']};
[fmdl,mat_idx]= ng_mk_cyl_models([1,1,.05],el_pos,[0.05,0,0.05], extra); 
show_fem(fmdl);
fmdl.stimulation(1).stim_pattern = [1;-1];
fmdl.stimulation(1).meas_pattern = [1,-1]; % dummy
img = mk_image(fmdl,1);
img.elem_data(mat_idx{2}) = 1e2;
img.fwd_solve.get_all_meas = 1;
vh=fwd_solve(img);

subplot(211)
show_fem(img);view(0,0);

return %%%%%%%%%%%%
% print_convert may14-resistive-target.png '-density 300 -resize 33%%'

imgv = rmfield(img,'elem_data');
imgv.node_data = vh.volt(:,1);

%imgv.calc_colours.clim = 0.10; % Colour limits
colours = calc_colours(imgv,[]);
subplot(212)
show_fem(fmdl); % STUPID MATlAB
hh= patch('Faces',fmdl.boundary,'Vertices',fmdl.nodes, 'facecolor','interp', ...
                'facevertexcdata',colours,'CDataMapping','direct'); 
view(0,0);
% print_convert may14-resistive-voltage.png '-density 300 -resize 33%%'

