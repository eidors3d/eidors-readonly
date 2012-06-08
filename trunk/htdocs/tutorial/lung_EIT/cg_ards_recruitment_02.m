%% SHOW EXAMPLE RESULTS FROM ONE STEP
load SUBJECT_1-c1.mat
EITCalcTimeSignal(db_c1.eitdata);
EITCalcFrequencySpectrum(db_c1.eitdata);
EITDisplayImages(db_c1.eitimages);
EITDisplayImages(db_c1.eitimages,1,'mtd');
% :clear
