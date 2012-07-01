close all
load SUBJECT_1-b.mat
load SUBJECT_1-c1.mat
load SUBJECT_1-c2.mat
load SUBJECT_1-c3.mat
load SUBJECT_1-c4.mat
load SUBJECT_1-d1.mat
load SUBJECT_1-d2.mat
load SUBJECT_1-d3.mat
load SUBJECT_1-d4.mat
clear matname
vars = who('db_*');
aggLungROI = eval([vars{1} '.eitimages.lungROI;']);
for i = 2:length(vars)
    aggLungROI = aggLungROI | eval([vars{i} '.eitimages.lungROI;']);
end
for i = 1:length(vars)
    eval([vars{i} '.eitimages.lungROI = aggLungROI;']);
    eval([vars{i} ...
        '.eitimages = EITCalcComplianceImage(' vars{i} '.eitimages);']);
end
clear i
limb_index = {'b','c1','c2','c3','c4'};

eitimages = eval(['db_' limb_index{1} '.eitimages']);
image_mask = eitimages.image_mask;
lungROI = eitimages.lungROI;
lindex = length(limb_index);
[row,col] = size(image_mask);
C = zeros(row,col,lindex);
pressure = zeros(lindex,1);
for i = 1:lindex
    C(:,:,i) = eval(['db_' limb_index{i} '.eitimages.complianceimage']);
    pressure(i) = eval(['db_' limb_index{i} '.eitimages.peep;']);
end
Ctab = TabulateImageData(C,image_mask);
lungROItab = TabulateImageData(lungROI,image_mask);
Ctab(lungROItab==0,:) = [];

Cmax = zeros(row,col);
Ctab = TabulateImageData(C,image_mask);
[Cmaxtab,index] = TabulateImageData(Cmax,image_mask);
Pstartab = Cmaxtab;
for i = 1:size(Ctab,1)
    if isinf(Ctab(i,1))
        Cmaxtab(i) = -Inf;
        Pstartab(i) = -Inf;
    else
        tempc = Ctab(i,:);
        cmax = max(tempc);
        pstaridx = find(tempc==cmax,1);
        pstar = eval(['db_' limb_index{pstaridx} '.eitimages.peep']);
        Cmaxtab(i)=cmax;
        Pstartab(i)=pstar;
    end
end
Cmax = ImageTableData(Cmaxtab,[row col],index);
Cmax(~image_mask) = -Inf;
Pstar = ImageTableData(Pstartab,[row col],index);
Pstar(~image_mask) = -Inf;
eitimages.Cmax  = Cmax;
eitimages.Pstar = Pstar;
for i = 1:lindex
    eval(['db_' limb_index{i} '.eitimages = EITCalcLungState(db_' ...
        limb_index{i} '.eitimages,Cmax,Pstar);']);
end
