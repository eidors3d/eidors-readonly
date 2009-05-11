function make_figures
   pig_data;



function [eelv,eilv] = pig_data
   if ~exist('p1130107.get','file');
      !wget http://eidors3d.sf.net/data_contrib/if-peep-acute-lung-injury/if_data_2003.zip
      !unzip if_data_2003.zip p1130107.get
   end
   vv= eidors_readdata('p1130107.get');
   eelv.p0  = vv(:,72);  eilv.p0  = vv(:,80);
   eelv.p20 = vv(:,554); eilv.p20 = vv(:,561);
   eelv.p30 = vv(:,795); eilv.p30 = vv(:,803);

