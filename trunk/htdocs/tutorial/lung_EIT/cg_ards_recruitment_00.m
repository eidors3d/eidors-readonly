zipfilecontents = unzip(...
    '../../data_contrib/cg-2012-ards-recruitment/cg_data_2012_p1.zip');
zipfilecontents = [zipfilecontents, unzip(...
    '../../data_contrib/cg-2012-ards-recruitment/cg_2012_ards_recruitment_code.zip')];
zipfilecontents{end+1} = 'zipfilecontents.mat';
save zipfilecontents.mat zipfilecontents
clear
