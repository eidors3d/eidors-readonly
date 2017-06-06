%delete('if_data_2003.zip');
for i = 1:length(zipfilecontents)
    delete(zipfilecontents{i});
end
