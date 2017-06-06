%delete('RPI_2d_chest_phantom_data.zip');
for i = 1:length(zipfilecontents)
    delete(zipfilecontents{i});
end
