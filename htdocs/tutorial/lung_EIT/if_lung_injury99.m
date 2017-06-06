%delete('if-experimental-lung-injury.zip');
for i = 1:length(zipfilecontents)
    delete(zipfilecontents{i});
end
