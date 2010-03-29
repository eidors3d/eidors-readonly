% $Id$

% Solve voltage for 3 different models
for idx=1:3
  if     idx==1; mdltype= 'd2C2';
  elseif idx==2; mdltype= 'd2d3c';
  elseif idx==3; mdltype= 'd2T2';
  end

  pat = 4; % Stimulation pattern to show

  imdl= mk_common_model(mdltype,16);
  img{idx} = calc_jacobian_bkgnd(imdl); 
  stim = mk_stim_patterns(16,1,'{ad}','{mono}',{'meas_current','rotate_meas'},-1);
  img{idx}.fwd_model.stimulation = stim(pat);
  img{idx}.fwd_solve.get_all_meas = 1;
end

% Show raw voltage pattern
for idx=1:3
  vh = fwd_solve(img{idx});
  imgn = rmfield(img{idx},'elem_data');
  imgn.node_data= vh.volt;
  subplot(2,3,idx); show_fem(imgn);
end
print -r125 -dpng backproj_solve02a.png;

% Calculate Equipotential lines
for idx=1:3
  vh = fwd_solve(img{idx});
  imgn = rmfield(img{idx},'elem_data');

  imgn.node_data= zeros(size(vh.volt,1),1);
  for v = 2:16
     imgn.node_data(vh.volt > vh.meas(v)) = v;
  end

  subplot(2,3,idx); show_fem(imgn);
end
print -r125 -dpng backproj_solve02b.png;
