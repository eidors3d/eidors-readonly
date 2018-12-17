fmdl= ng_mk_ellip_models([4,0.8,1.1,.5],[16,1.7,2.3],[0.05]);
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip4{:});

% "Square" electrode layout
idx = reshape(1:32,2,[])';
idx(2:2:end,:) = fliplr(idx(2:2:end,:));
extraflip= [4:12]; % This belt was made slightly differently
idx(extraflip,:) = fliplr(idx(extraflip,:));
fmdl.electrode(idx) = fmdl.electrode(:);
