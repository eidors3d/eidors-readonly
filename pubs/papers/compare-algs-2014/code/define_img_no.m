function [iname, Dt] =  define_img_no(i,expis,insps)
 
   if nargin==1;
      expis= zeros(0,6);
      insps= zeros(0,6);
   end

   switch i;      
     case  1; iname = 'VT,EIT,ZEEP,21,#1';
              Dt= [expis(:,1),insps(:,1)];
     case  2; iname = 'VT,EIT,PEEP,21,#1';
              Dt= [expis(:,2),insps(:,2)];
     case  3; iname = 'VT,EIT,ZEEP,100';
              Dt= [expis(:,3),insps(:,3)];
     case  4; iname = 'VT,EIT,PEEP,100';
              Dt= [expis(:,4),insps(:,4)];
     case  5; iname = 'VT,EIT,ZEEP,21,#2';
              Dt= [expis(:,5),insps(:,5)];
     case  6; iname = 'VT,EIT,PEEP,21,#2';
              Dt= [expis(:,6),insps(:,6)];
     case 10; iname = '???';
              Dt= [expis(:,2),insps(:,2),expis(:,1),insps(:,1)];
     case 11; iname = '???';
              Dt= [expis(:,4),insps(:,4),expis(:,3),insps(:,3)];
     case 12; iname = 'CEELV,EIT,Z-P,21,#1';
              Dt= expis(:,[1,2]);
     case 13; iname = 'CEELV,EIT,Z-P,21,#2';
              Dt= expis(:,[5,6]);
     case 14; iname = 'CEELV,EIT,Z-P,100';
              Dt= expis(:,[3,4]);
     case 15; iname = 'CEELV,EIT,Z-Z,21';
              Dt= expis(:,[1,5]);
     case 16; iname = 'CEELV,EIT,P-P,21';
              Dt= expis(:,[2,6]);
     case 17; iname = 'CEELV,EIT,Z21-Z100,#1';
              Dt= expis(:,[1,3]);
     case 18; iname = 'CEELV,EIT,Z21-Z100,#2';
              Dt= expis(:,[5,3]);
     case 19; iname = 'CEELV,EIT,P21-P100,#1';
              Dt= expis(:,[2,4]);
     case 20; iname = 'CEELV,EIT,P21-P100,#2';
              Dt= expis(:,[6,4]);
     otherwise; error('value of i not recognized (%d)',i);
   end

