function param= aa_resetup(reponse,r1,r2,r3)
% RESETUP: prepare fwd_model for EIT image reconstruction
% prepare les variables necessaire pour 
% fonctionner res.m, grille.m, et imgr.m dynproj.m
%
% format  res_setup(' (abcd)(012345)[-|+][x] ',
%                     position_initial,
%                    ratio elliptique [x/y] );
%
% if string begins with - then no verbose output
%
% position_initial = [0 2] DEFAUT
% [x] --> ne pas calculer DVV
%
% (c) 1993-2002 A.Adler
% $Id: aa_resetup.m,v 1.1 2004-07-20 00:39:50 aadler Exp $

VERBOSE=1;

if nargin==0
  reponse=input('quelle geometrie,position, etc.(''?'' pour aide): ','s'); 

  if reponse=='?'
    fprintf(['preparer quelle geometrie\n' ...
             ' a: tank with 120 elements    0: circ tank.\n' ...
             ' b: tank with 216 elements    1: thorax lev1 \n' ...
             ' c: tank with 64 elements     2: thorax lev2 \n' ...
             ' d: tank with 256 elements    3: thorax lev3 \n' ...
             ' e: tank with 576 elements    4: thorax lev4 \n' ...
             ' f: tank with 1024 elements   5: thorax lev5 \n' ...
             ' i: tank with 352 elements    9: thorax of dog\n'... 
             ' j: tank with 86 elements     \n'... 
             ' m: bassin 3D a 6*64 elem     \n' ...
             ' n: bassin 3D a 6*256 elem    \n' ...
             ' o: bassin 3D a 6*576 elem    \n' ...
             ' p: bassin 3D a 6*1024 elem   \n' ...
             ' s: bassin 3D a 12*64 elem    \n' ...
             ' t: bassin 3D a 12*256 elem   \n' ...
             '\nexemple> a0,0 1 16,1' ...
             ' a0= geometrie,0 1 16=[sce puit #electr],1= elliptique' ]); 
    reponse=input('votre choix: ','s'); 
  end

elseif nargin==2
  reponse=[reponse ',' r1];
elseif nargin==3
  reponse=[reponse ',' r1 ',' r2];
elseif nargin==4
  reponse=[reponse ',' r1 ',' r2 ',' r3];
end

if reponse(1)=='-'
  VERBOSE=0;
  reponse = reponse(2:length(reponse));
end

i=[ find(reponse==','), find(reponse=='/')];

if length(i)==0
  ellip=[];
  pos_i=[];
elseif length(i)==1
  ellip=[];
  pos_i=sscanf(reponse(i+1:length(reponse)),'%f')';
  reponse=reponse(1:i-1);
elseif length(i)==2
  ellip=sscanf(reponse(i(2)+1:length(reponse)),'%f')';
  pos_i=sscanf(reponse(i(1)+1:i(2)-1),'%f')';
  reponse=reponse(1:i(1)-1);
else
  error('parametres non comprises');
end


  if isempty(ellip)
    if reponse(2)=='9'
      ellip=1.2;
    else
      ellip=1;
    end;
  end
  if VERBOSE ; fprintf('ellip= %2.2f\n',ellip); end

  if length(pos_i)==0
    if any( reponse(2)== '123459' )
      pos_i= [8 9];
    else
      pos_i= [0 1];
    end %if reponse
    elec=16;
  elseif length(pos_i)==1
    pos_i= pos_i+[0 1];
    elec=16;
  elseif length(pos_i)==2
    elec= 16;
  elseif length(pos_i)==3
    elec= pos_i(3);
    if ~any(elec==[4 8 16 32 64])
      error('Nombre d''electrodes invalide');
    end %if ~any
    pos_i= pos_i(1:2);
  else
    error('POS_I non interpretable');
  end %if length

  if VERBOSE;
  fprintf('[sce puit #electr] = [%2.1f  %2.1f  %2d]\n',pos_i(1),pos_i(2),elec);
  end
  if exist('OCTAVE_VERSION')
  choix=[reponse(1:2), ...
         setstr( pos_i+toascii('0')*(pos_i<10)+(toascii('a')-10)*(pos_i>9) )];
  else
  choix=[reponse(1:2), ...
         setstr( pos_i+abs('0')*(pos_i<10)+(abs('a')-10)*(pos_i>9) )];
  end
  if strcmp(ChoiX,choix)
    if VERBOSE ;
        disp([ChoiX '  ' choix])
        disp(['Resetup deja fait pour choix= ' ChoiX]);
    end
    return;
  else
    ChoiX= choix;
    if VERBOSE ; disp(['choix= ' ChoiX]); end
    clear global ZparaM ZproJ
  end

  if reponse(1)=='a' | reponse(1)=='b';
    phi=(2*pi/64)*(0:63);
    NODE=[ [0;0] ...
      .2*[ sin(phi(1:16:64)); cos(phi(1:16:64)) ] ...
      .4*[ sin(phi(1:8:64)); cos(phi(1:8:64)) ] ... 
      .6*[ sin(phi(1:4:64)); cos(phi(1:4:64)) ] ...
      .8*[ sin(phi(3:4:64)); cos(phi(3:4:64)) ] ... 
     .95*[ sin(phi(1:2:64)); cos(phi(1:2:64)) ] ];
    ELEM= [ [2:5; 3:5 2; ones(1,4)] ...
            [2:5; 3:5 2; 7:2:13 ] ...
            [6:13; 7:13 6; 2+ceil(0:.5:3) 2] ...      
            [6:13; 7:13 6; 15:2:29 ] ...
            [14:29;15:29 14;6+ceil(0:.5:7) 6] ...
            [14:29;15:29 14;30:45] ... 
            [45 30:44;30:45;14:29] ...
            [45 30:44;30:45;46:2:77] ...
            [46:77;47:77 46;30+floor(0:.5:15.5)] ]; 
    MES=[46:32/elec:77];     
    CPTR=[1:72 73*ones(1,48)];

    if reponse(1)=='b'
      NODE=[ NODE [sin(phi); cos(phi)] ];   
      ELEM=[ ELEM ...
           [46:77;47:77 46;79:2:142] ...
           [78:141;79:141 78;46+ceil(0:.5:31) 46] ];
      MES= [78:64/elec:141];
      CPTR= [1:120 121*ones(1,96)];
    end %if reponse=='b'
  end %if reponse=='a' | reponse=='b';

  if reponse(1)=='j' 
    thomas_circ
  end

  if any(reponse(1)=='mnop')
     niveaux= [-.1;0;.1];
     reponse(1)= reponse(1)-abs('m')+abs('c');
     if reponse(2) ~= '0'
       niveaux=niveaux*100;
     end
  end %if reponse='mnop'

  if any(reponse(1)=='st')
     niveaux= [-.25;-.1;0;.1;.25];
     reponse(1)= reponse(1)-abs('s')+abs('c');
     if reponse(2) ~= '0'
       niveaux=niveaux*100;
     end
  end %if reponse='st'

  if any(reponse(1)=='cdefi')
    if exist('OCTAVE_VERSION')
        N= ( toascii(reponse(1))-toascii('b') )*4;
    else
        N= ( abs(reponse(1))-abs('b') )*4;
    end
    if reponse(1)=='i', N=8; end;
    ELEM=[];
    NODE= [0;0];
    int=1;
    for k=1:N
      phi= (0:4*k-1)*pi/2/k;
      NODE= [NODE k/N*[sin(phi);cos(phi)]];

      ext= 2*(k*k-k+1);
      idxe=[0:k-1; 1:k];
      idxi=[0:k-1]; 
      elem= [ ext+idxe ext+2*k+[-idxe idxe] ...
                      ext+rem(4*k-idxe,4*k) ...
              ext+idxe ext+2*k+[-idxe idxe] ...
                      ext+rem(4*k-idxe,4*k);
              int+idxi int+2*(k-1)+[-idxi idxi] ... 
                          int+rem(4*(k-1)-idxi, ...
                               4*(k-1)+(k==1) ) ...
              ext+4*k+1+idxi ext+6*k+ ...
               [1-idxi 3+idxi] ext+8*k+3-idxi ];
      for j=1:k
        r1= rem(j+k-1,3)+1;
        r2= rem(j+k,3)+1;
        r3= 6-r1-r2;
        elem([r1 r2 r3],j+k*(0:7) )= elem(:,j+k*(0:7));
      end
  
      ELEM=[ ELEM elem(:,1:(8-4*(k==N))*k) ];
      int=ext;
    end %for k=1:N
  
    MES= ext+N*4/elec*([0:elec-1]);
 
    CPTR= [1:(2*N-2)^2 ((2*N-2)^2+1)*ones(1,8*N-4) ];
  end %if reponse= 'c' 'd' 'e'

  if reponse(1)=='i'
    phi=(2*pi/64)*(0:63);
    NODE=[ .95*NODE [sin(phi); cos(phi)] ];
    ELEM=[ ELEM ...
       [114:145;115:145 114;147:2:209] ...
       [146:209;147:209 146;114+ceil(0:.5:31) 114] ];
    MES= [146:64/elec:209];
    CPTR= [1:256 257*ones(1,96)];
  end %if reponse=='i'

  j=pos_i(1); k=floor(j);
  MES=[MES( 1+rem( (0:elec-1)+k,elec) )+floor(MES(1:2)*[-1;1]*(j-k));
       MES( 1+rem( (1:elec)+k,elec) )+floor(MES(1:2)*[-1;1]*(j-k)) ];

  cour= [ MES(1,:)' MES(1,1+rem(elec+(0:elec-1)+pos_i*[-1;1],elec) )'; 
           ones(elec,1)*[-1 1] ];

  volt= [1; zeros(elec,1)];

  ELS=rem(rem(0:elec^2-1,elec)-floor((0:elec^2-1)/elec)+elec,elec)';
  ELS=~any(rem( elec+[-1 0 [-1 0]+pos_i*[-1;1] ] ,elec)' ...
		     *ones(1,elec^2)==ones(4,1)*ELS')';

  if any(reponse(2)=='12345')
    niv= reponse(2)-abs('0');
    x_coord= [ ...
      0,60,123,173,223,202,144,75,0,-75,-144,-202,-223,-173,-123,-60;
      0,50,105,138,144,144,109,50,0,-50,-109,-144,-144,-138,-105,-50;
      0,52, 99,133,148,141,110,61,0,-61,-110,-141,-148,-133,- 99,-52;
      0,51, 92,129,148,136, 96,47,0,-47,- 96,-136,-148,-129,- 92,-51;
      0,49, 92,128,148,141,111,64,0,-64,-111,-141,-148,-128,- 92,-49 ];
    y_coord= [ ...
      123,116, 91,42,-4,-67,-105,-119,-108,-119,-105,-67,-4,42, 91,116;
      129,132,112,62, 3,-57,-101,-110,-107,-110,-101,-57, 3,62,112,132;
%     128,132,112,62, 3,-57,-101,-110,-102,-110,-101,-57, 3,62,112,132;
      116,112, 92,53, 3,-48,- 88,-106,-105,-106,- 88,-48, 3,53, 92,112;
      143,130, 99,63,14,-35,- 68,- 82,- 82,- 82,- 68,-35,14,63, 99,130;
      136,128,103,68,23,-25,- 62,- 78,- 80,- 78,- 62,-25,23,68,103,128 ];
    geo= [x_coord(niv,[13:16 1:12])',  ...
              y_coord(niv,[13:16 1:12])'];
    a_max= size(geo,1);
    ab_geo=sqrt(sum(([ geo; geo(1,:) ]').^2)');
    nn= zeros(size(NODE));
    for i=1:size(NODE,2);
      angle = rem(a_max*atan2( NODE(2,i), ...
            NODE(1,i) )/2/pi+a_max,a_max)+1;
      fac=(  (floor(angle+1.001)- angle)* ...
              ab_geo(floor(angle+.001)) + ...
             (angle-floor(angle+.001))* ...
              ab_geo(floor(angle+1.001))  );
      nn(:,i)= NODE(:,i)* fac;
    end  %for i=1:size
    NODE=nn;
  end %if reponse(2)

  NODE(1,:)= NODE(1,:)*sqrt(2)/sqrt(1+ellip^2);
  NODE(2,:)= NODE(2,:)*sqrt(2)*ellip/sqrt(1+ellip^2);


d=  size(ELEM,1);       %dimentions+1
n= size(NODE,2);        %NODEs
e= size(ELEM,2);        %ELEMents     

if length(reponse)>2 if reponse(3)=='+'
  disp('CPTR adapte'); end
else
  CPTR= 1:e;
  if VERBOSE; disp('CPTR = 1:elements'); end
end; 
 

fichier= ['eptr_' reponse(1:2) '.mat' ];
if exist(fichier)==2 

  load(fichier)
  disp([ fichier ' charge']);

else %if exist

 if VERBOSE; fprintf('creant EPTR ..'); end
 tic=cputime;
 npx= 64; npy= 64;    %nombre de points en chaque direction
 %used 
 [x y]=meshgrid( ...
     linspace( min(NODE(1,:))*1.05, max(NODE(1,:))*1.05 ,npx ), ...
    -linspace( min(NODE(2,:))*1.05, max(NODE(2,:))*1.05 ,npy )  ); 
 v_yx= [-y(:) x(:)];
 tourne= [0 -1 1;1 0 -1;-1 1 0];
 EPTR=zeros(npy,npx);
 for j= 1:e
   xy= NODE(:,ELEM(:,j))';
   a= xy([2;3;1],1).*xy([3;1;2],2)- xy([3;1;2],1).*xy([2;3;1],2);
   endr=find( y(:)<=max(xy(:,2)) & y(:)>=min(xy(:,2)) ...
	    & x(:)<=max(xy(:,1)) & x(:)>=min(xy(:,1)) );
   aa= sum(abs(ones(length(endr),1)*a'+ ...
	       v_yx(endr,:)*xy'*tourne)');
   endr( abs( (abs(sum(a))-aa) ./ sum(a)) >1e-8)=[];
   EPTR(endr)= j;
 end %for j=1:ELEM
 if VERBOSE ; disp(cputime-tic); end

end % if ~exist

if exist('niveaux')
  node=NODE;
  elem= [ELEM([1 1 2 3 ],:) ...
         ELEM([2 1 2 3],:) ELEM([3 2 1 3],:)]; 
  NODE= [node; niveaux(1)*ones(1,n) ];
  ELEM= [];
 
  for k=2:length(niveaux);
    NODE=[NODE  [node; niveaux(k)*ones(1,n)] ];
    ELEM= [ELEM (elem + ...
       [[(k-1)*n*ones(1,e);(k-2)*n*ones(3,e)] ...
        [(k-1)*n*ones(2,e);(k-2)*n*ones(2,e)] ...  
        [(k-1)*n*ones(3,e);(k-2)*n*ones(1,e)]] ) ];
  end %for k

  MES= MES + floor(length(niveaux)/2)*n;
  CPTR= [CPTR CPTR CPTR CPTR CPTR CPTR];
  cour(1:elec,:)= cour(1:elec,:)+ floor(length(niveaux)/2)*n;

end %if exist('niveaux')

d=  size(ELEM,1);       %dimentions+1
n= size(NODE,2);        %NODEs
e= size(ELEM,2);        %ELEMents     
p= size(volt,1)-1;
c= max(CPTR);

AIRE=zeros(e,1);
for i=1:e
  AIRE(i)= abs(det([ones(1,d);NODE(:,ELEM(:,i))]))/(d-1)/(d-2);
end

CC= sparse((1:d*e),ELEM(:),ones(d*e,1), d*e, n);

if 0;
    SS= spalloc( d*e , d*e, d*d*e);
    for j=1:e
      a=  inv([ ones(d,1) NODE( :, ELEM(:,j) )' ]);
      SS(d*(j-1)+1:d*j,d*(j-1)+1:d*j)=  ...
           2*a(2:d,:)'*a(2:d,:)/(d-1)/(d-2)/abs(det(a));
    end %for j=1:ELEMs 
else
    SSiidx= floor([0:d*e-1]'/d)*d*ones(1,d) + ones(d*e,1)*(1:d) ;
    SSjidx= [1:d*e]'*ones(1,d);
    SSdata= zeros(d*e,d);
    dfact= (d-1)*(d-2); % Note this wont work for d>3
    for j=1:e
      a=  inv([ ones(d,1), NODE( :, ELEM(:,j) )' ]);
      SSdata(d*(j-1)+1:d*j,1:d)= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));
    end %for j=1:ELEMs 
    SS= sparse(SSiidx,SSjidx,SSdata);
end

QQ=sparse(cour(1:p,:),(1:p)'*ones(1,size(cour,2)), ...
                     cour(p+1:2*p,:),n,p );

if any(reponse(2:length(reponse))=='x')
    disp('je ne calcule pas DVV');
else
  if VERBOSE ; fprintf('creant DVV ..'); end
  tic=cputime;
  i=1:n; i(volt(1,:))=[];
  DVV= zeros(sum(ELS),c);
  z= CC'*SS*CC;
  zinv=zeros(n);
  if exist('OCTAVE_VERSION')
      zinv(i,i)= spinv(z(i,i));
  else
      zinv(i,i)= inv(z(i,i));
  end
  dzi= zinv(MES(1,:),:)- zinv(MES(2,:),:);
  vvh= dzi*QQ;
  sv= zinv*QQ;
  for j=1:c
    idx= find(ones(d,1)*(CPTR==j));
    dq= dzi*(  CC(idx,:)'*SS(idx,idx)*CC(idx,:)  )*sv;
    DVV(:,j)= dq(ELS)./vvh(ELS);
  end
  if VERBOSE ; disp(cputime-tic); end
end

   
if 1 %nargin==0
% grille; %(0);

  %s=.5* (0:1/32:1).^(.5);
  %s=[ .5-fliplr(s(2:33)) s+.5 1]'*ones(1,3);

  %s=[ zeros(1,10) (0:1/44:1) ones(1,11) ]'*ones(1,3);
  %s=[ ones(1,10) (1:-1/44:0) zeros(1,10) 1 ]'*ones(1,3);

  s= [flipud(fliplr(hot(64))) ;hot(64)];
  s= [s([1:2:63 64:2:128],:)*.7+.3;[.5 .3 .6]];
  colormap(s);
end %if nargin

param.nodes = NODE;
param.elems = ELEM;
param.MES   = MES;
param.QQ    = QQ;
param.CC    = CC;
param.SS    = SS;
param.EPTR  = EPTR;
param.ELS   = ELS;
param.CPTR  = CPTR;
param.DVV   = DVV;
param.AIRE  = AIRE;
