% PROCESS EACH STEP OF THE PROTOCOL
basename = './DATA';
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
