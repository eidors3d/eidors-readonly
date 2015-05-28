data.imageRate = 13;

data.patient.ROI.Inside =thorax_ROI*100; % to scale it up to 100
data.patient.ROI.RightLung =rlung_ROI*100;
data.patient.ROI.LeftLung =llung_ROI*100;
data.patient.ROI.Heart =zeros(size(imgs,1),size(imgs,2));

% put to dummy because they are missing
data.patient.halfChest = 'NaN';
data.patient.height = 'NaN';
data.patient.weight = 'NaN';
data.patient.gender = 'NaN';

data.measurement.Position.transversal = zeros (1,size(imgs,3));
data.measurement.Position.longitudinal = zeros (1,size(imgs,3));
data.measurement.ImageQuality = 100*ones(1,size(imgs,3));
data.measurement.ElectrodeQuality = zeros(size(imgs,3),32);
data.measurement.ZeroRef = imgs;

data.injctionPattern= 'NaN';
data.SensorBelt.NumEl= 'NaN';

data.measurement.CompositValue=squeeze(sum(sum(imgs,2),1));

save('file-for-IBEX.mat','data');

