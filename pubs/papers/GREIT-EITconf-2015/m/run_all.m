d = dir('../data/Skip*');
d = d([d.isdir]);
for i = 1:numel(d);
   tank_recon(d(i).name);
end