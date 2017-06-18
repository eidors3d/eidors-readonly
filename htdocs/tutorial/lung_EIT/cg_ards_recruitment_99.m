load zipfilecontents.mat
delete(zipfilecontents{:});
delete('SUBJECT*.mat');
rmdir('DATA','s'); % remove the data dir
