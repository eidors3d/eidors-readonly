function process_and_Save(basename, filename, range, maneuver, PEEP, dP)
   close all
   [eitdata,s] = EITReadData(basename, filename);
   eitdata.maneuver = maneuver;
   eitdata.peep = PEEP;
   eitdata.deltap = dP;
   if isempty(range)
       range = [1,size(eitdata.voltage_data,2)];
   end
   [eitdata,s]=EITFilterData(eitdata,'bandpass',range);
   [eitimages,s]=EITReconstructImages(eitdata);
   [eitdata,eitimages,s]=EITCalcTidalImages(eitdata,eitimages,range);
   
   % Temp line for counting breaths
   eitdata.tidalindices = eitdata.tidalindices;
   [eitimages,p,s]=EITCalcLungRoi(eitimages);
   [eitimages,s]=EITCalcComplianceImage(eitimages);
   
   [tok,rem]=strtok(eitimages.file_name,'_');
   tok = strtok(rem,'.');
   tok(1)=[];
   matname = [eitimages.subject_id '-' tok '.mat']
   f1='eitimages';f2='eitdata';f3='filename';f4='matname';
   eval(['db_' tok ...
       ' = struct(f1,eitimages,f2,eitdata,f3,filename,f4,matname);']);
   
   clear PEEP dP maneuver range s ans tok rem f1 f2 f3 f4 p
   clear eitdata eitimages filename basename

   save(matname)
