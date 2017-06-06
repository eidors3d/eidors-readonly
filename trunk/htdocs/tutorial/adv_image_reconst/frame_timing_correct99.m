%delete('if-neonate-ventilated.zip');
for i = 1:length(zipfilecontents)
    delete(zipfilecontents{i});
end
