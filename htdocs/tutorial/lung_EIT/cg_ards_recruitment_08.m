close all
for i = 1:lindex
    EITDisplayImages(eval(['db_' limb_index{lindex-i+1} '.eitimages']),1,'collapseOD');
end
