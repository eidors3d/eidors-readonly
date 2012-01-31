% ----------------- TUTORIAL: REGIONAL LUNG MECHANICS ---------------------
%
% AUTHORS: Camille Gomez-Laberge, John H Arnold & Gerhard K Wolf
%
% DATASET DESCRIPTION: Stepwise recruitment maneuver and PEEP titration of
% a patient with acute respiratory distress syndrome. Performed at
% Children's Hosptial Boston in 2009.
%
% PLEASE CITE THE FOLLOWING ARTICLE WHEN USING THIS SOFTWARE:
%   Gomez-Laberge C, Arnold JH & Wolf GK. A unified approach for EIT
%   imaging of regional overdistension and atelectasis in acute lung
%   injury. IEEE Trans Med Imag, VOL: P-P, 2012.
%
% SOFTWARE CONTRIBUTED TO EIDORS UNDER THE FOLLOWING PROJECTS
%   RegionalLungMechanics/
% WITH MODIFICATIONS OF THE FOLLOWING GRAPHICS FUNCTIONS
%   eidors/graphics_matlab/calc_colours.m
%   eidors/graphics_matlab/eidors_colours.m
%   eidors/graphics_matlab/scale_for_display.m
%   Note: modified sections begin and end with the tag %|||RLM edit
%
% -------------------------------------------------------------------------
clear
%% ------------------------------ STAGE 1 ---------------------------------
% PROCESS EACH STEP OF THE PROTOCOL
basename = 'E:/TUTORIAL/DATA';
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_b.get'; 
range =[]; maneuver='increment'; PEEP=14; dP=5;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_c1.get'; 
range =[]; maneuver='increment'; PEEP=15; dP=15;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_c2.get'; 
range =[]; maneuver='increment'; PEEP=20; dP=15;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_c3.get'; 
range =[]; maneuver='increment'; PEEP=25; dP=15;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_c4.get'; 
range =[]; maneuver='increment'; PEEP=30; dP=15;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_d1.get'; 
range =[1,765]; maneuver='decrement'; PEEP=20; dP=6;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_d2.get'; 
range =[]; maneuver='decrement'; PEEP=18; dP=6;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_d3.get'; 
range =[368,780]; maneuver='decrement';  PEEP=16; dP=5;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
filename = 'STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_d4.get'; 
range =[411,780]; maneuver='decrement'; PEEP=14; dP=4;
process_and_save(basename, filename, range, maneuver, PEEP, dP);
    
%% SHOW EXAMPLE RESULTS FROM ONE STEP
load SUBJECT_1-c1.mat
EITCalcTimeSignal(db_c1.eitdata);
EITCalcFrequencySpectrum(db_c1.eitdata);
EITDisplayImages(db_c1.eitimages);
EITDisplayImages(db_c1.eitimages,1,'mtd');
clear
% -------------------------------------------------------------------------

%% ------------------------------ STAGE 2 ---------------------------------
% AGGREGATE ALL STEPS IN THE RECRUITMENT MANEUVER
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
vars = who;
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

%% SHOW MAX COMPLIANCE AND PSTAR MAPS
EITDisplayImages(eitimages,1,'Cmax')
EITDisplayImages(eitimages,1,'Pstar')

%% SHOW OVERDISTENSION AND ATELECTASIS MAPS
close all
for i = 1:lindex
    EITDisplayImages(eval(['db_' limb_index{lindex-i+1} '.eitimages']),1,'collapseOD');
end

%% ------------------------------ STAGE 3 ---------------------------------
% AGGREGATE ALL STEPS IN THE PEEP TITRATION MANEUVER
close all
clear
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
vars = who;
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

limb_index = {'d1','d2','d3','d4'};

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

%% SHOW MAX COMPLIANCE AND PSTAR MAPS
EITDisplayImages(eitimages,1,'Cmax')
EITDisplayImages(eitimages,1,'Pstar')

%% SHOW OVERDISTENSION AND ATELECTASIS MAPS
close all
for i = 1:lindex
    EITDisplayImages(eval(['db_' limb_index{lindex-i+1} '.eitimages']),1,'collapseOD');
end
