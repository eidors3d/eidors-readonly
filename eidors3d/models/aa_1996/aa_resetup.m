function param= aa_resetup(varargin)
% RESETUP: prepare fwd_model for EIT image reconstruction
%
% format  res_setup(' (abcd)(012345)[-|+][x] ',
%                     position_initial,
%                    ratio elliptique [x/y] );
%
% if string begins with - then no verbose output
%
% position_initial = [0 2] DEFAUT
%
% (c) 1993-2002 A.Adler
% $Id: aa_resetup.m,v 1.4 2004-07-24 04:11:07 aadler Exp $

inp = proc_input_args( varargin{:} );

niveaux= [];
if inp.selec >= 'c' && inp.selec <= 'f'
 % 2D models
   rings    = (inp.selec-abs('c')+1 )*4;
elseif  inp.selec=='i'
   rings    =  8;
elseif  inp.selec>='k'
  % 3D Models

  if inp.selec >= 'k' && inp.selec <= 'n'
     niveaux= [-.1;0;.1];
     rings    = (inp.selec-abs('k')+1)*4;
  elseif inp.selec >= 'o' && inp.selec <= 'r'
     niveaux= [-.25;-.1;0;.1;.25];
     rings    = (inp.selec-abs('o')+4)*4;
  end 

  % if using thorax model, then we need to expand size of layers
  if inp.reponse(2) ~= '0'
    niveaux=niveaux*100;
  end

end %if inp.selec > 'c'


if inp.selec=='a' || inp.selec=='b';
  mdl= mk_2d_fem_type1( inp.selec, inp.n_elec );
else
  mdl= mk_circ_tank(rings, niveaux, inp.n_elec, inp.n_planes );
end

%  if inp.selec=='i'
%    phi=(2*pi/64)*(0:63);
%    NODE=[ .95*NODE [sin(phi); cos(phi)] ];
%    ELEM=[ ELEM ...
%       [114:145;115:145 114;147:2:209] ...
%       [146:209;147:209 146;114+ceil(0:.5:31) 114] ];
%    MES= [146:64/elec:209];
%    CPTR= [1:256 257*ones(1,96)];
%  end %if reponse=='i'

param= mdl;
param.type= 'model';
return;

  j=pos_i(1); k=floor(j);
  MES=[MES( 1+rem( (0:elec-1)+k,elec) )+floor(MES(1:2)*[-1;1]*(j-k));
       MES( 1+rem( (1:elec)+k,elec) )+floor(MES(1:2)*[-1;1]*(j-k)) ];

  cour= [ MES(1,:)' MES(1,1+rem(elec+(0:elec-1)+pos_i*[-1;1],elec) )'; 
           ones(elec,1)*[-1 1] ];

  volt= [1; zeros(elec,1)];

  ELS=rem(rem(0:elec^2-1,elec)-floor((0:elec^2-1)/elec)+elec,elec)';
  ELS=~any(rem( elec+[-1 0 [-1 0]+pos_i*[-1;1] ] ,elec)' ...
		     *ones(1,elec^2)==ones(4,1)*ELS')';



d= size(ELEM,1);       %dimentions+1
n= size(NODE,2);        %NODEs
e= size(ELEM,2);        %ELEMents     


fichier= ['eptr_' reponse(1:2) '.mat' ];
if exist(fichier)==2 

  load(fichier)
  disp([ fichier ' charge']);

else %if exist

 eidors_msg('creating EPTR ..',5);
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
  cour(1:elec,:)= cour(1:elec,:)+ floor(length(niveaux)/2)*n;

end %if exist('niveaux')

d=  size(ELEM,1);       %dimentions+1
n= size(NODE,2);        %NODEs
e= size(ELEM,2);        %ELEMents     
p= size(volt,1)-1;

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
  eidors_msg('creating jacobian ..',4);
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
    idx= find(ones(d,1)*((1:e)==j));
    dq= dzi*(  CC(idx,:)'*SS(idx,idx)*CC(idx,:)  )*sv;
    DVV(:,j)= dq(ELS)./vvh(ELS);
  end
end


param.nodes = NODE;
param.elems = ELEM;
param.MES   = MES;
param.QQ    = QQ;
param.CC    = CC;
param.SS    = SS;
param.EPTR  = EPTR;
param.ELS   = ELS;
param.DVV   = DVV;
param.AIRE  = AIRE;

% process input args to select parameters to use
% Output is:
%   inparam.VERBOSE=     VERBOSE;
%   inparam.ellip=       ellip;
%   inparam.pos_i=       pos_i;
%   inparam.reponse=     reponse;
%   inparam.n_elec=      elec;
%   inparam.n_planes=    elec_planes;
%   inparam.selec=       reponse(1);
function inparam= proc_input_args( reponse, r1, r2, r3 )

if nargin==0
  reponse=input('Choose model geometry. (type "?" for help): ','s'); 

  if reponse=='?'
    fprintf(['Choose a model geometery to prepare\n' ...
             ' a: 2D tank with 120 elements  0: circ tank.\n' ...
             ' b: 2D tank with 216 elements  1: thorax lev1 \n' ...
             ' c: 2D tank with 64 elements   2: thorax lev2 \n' ...
             ' d: 2D tank with 256 elements  3: thorax lev3 \n' ...
             ' e: 2D tank with 576 elements  4: thorax lev4 \n' ...
             ' f: 2D tank with 1024 elements 5: thorax lev5 \n' ...
             ' i: 2D tank with 352 elements  9: thorax of dog\n'... 
             ' k: 3D tank with 6*64 elem     \n' ...
             ' l: 3D tank with 6*256 elem    \n' ...
             ' m: 3D tank with 6*576 elem    \n' ...
             ' n: 3D tank with 6*1024 elem   \n' ...
             ' o: 3D tank with 12*64 elem    \n' ...
             ' p: 3D tank with 12*256 elem   \n' ...
             ' q: 3D tank with 12*576 elem    \n' ...
             ' r: 3D tank with 12*1024 elem   \n' ...
             '\nexample> a0,0 1 16 1' ...
             ' a0= geometry,0 1 16=[source sink #electr #planes]' ]); 
    reponse=input('Your choice: ','s'); 
  end

elseif nargin==2
  reponse=[reponse ',' r1];
elseif nargin==3
  reponse=[reponse ',' r1 ',' r2];
elseif nargin==4
  reponse=[reponse ',' r1 ',' r2 ',' r3];
end

VERBOSE=1;

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
  error('Input parameters not understood');
end

  if isempty(ellip)
    if reponse(2)=='9'
      ellip=1.2;
    else
      ellip=1;
    end;
  end
  eidors_msg('Elliptical excentricity (1=circular)= %2.2f',ellip,4);

  if length(pos_i)>=4
    elec_planes= pos_i(4);
  else 
    elec_planes= 1;
  end

  if length(pos_i)==0
    pos_i= [0 1];
    elec=16;
  elseif length(pos_i)==1
    pos_i= pos_i+[0 1];
    elec=16;
  elseif length(pos_i)==2
    elec= 16;
  elseif length(pos_i)<=3
    elec= pos_i(3);
    if ~any(elec==[4 8 16 32 64])
      error('Invalid number of electrodes');
    end %if ~any
    pos_i= pos_i(1:2);
  end %if length

  eidors_msg('[source sink #electr #planes] = [%2d %2d %2d %2d]', ...
             pos_i(1),pos_i(2),elec, elec_planes, 1);

inparam.VERBOSE=     VERBOSE;
inparam.ellip=       ellip;
inparam.pos_i=       pos_i;
inparam.reponse=     reponse;
inparam.n_elec=      elec;
inparam.n_planes=    elec_planes;
inparam.selec=       reponse(1);

function param = mk_2d_fem_type1( selec, n_elec )
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
    MES=[46:32/n_elec:77];     
    bdy_nodes=[ 46:77;47:77,46 ];

    if selec=='b'
      NODE=[ NODE [sin(phi); cos(phi)] ];   
      ELEM=[ ELEM ...
           [46:77;47:77 46;79:2:142] ...
           [78:141;79:141 78;46+ceil(0:.5:31) 46] ];
      MES= [78:64/n_elec:141];
    end %if reponse=='b'
    bdy_nodes=[ 78:141;79:141,78 ];

    param.nodes = NODE';
    param.elems = ELEM';
    param.boundary = bdy_nodes';
    param.gnd_node = 1; % node at center of the tank
    param.name= sprintf('2D FEM for EIT with %d nodes',size(NODE,2));

% create point electrodes
    for i=1:length(MES)
       param.electrode(i).nodes= MES(i);
       param.electrode(i).z_contact= 0;
    end
   

% reshape model based on elliptical excentricity and a thorax level
function newmdl = reshape_mdl( mdl, niv, ellip );
% mdl is input model
% niv is the thorax level (1 to 5)
% ellip is elliptical excentricity
    NODE= mdl.nodes';
    ELEM= mdl.elems';
    x_coord= [ ...
      0,60,123,173,223,202,144,75,0,-75,-144,-202,-223,-173,-123,-60;
      0,50,105,138,144,144,109,50,0,-50,-109,-144,-144,-138,-105,-50;
      0,52, 99,133,148,141,110,61,0,-61,-110,-141,-148,-133,- 99,-52;
      0,51, 92,129,148,136, 96,47,0,-47,- 96,-136,-148,-129,- 92,-51;
      0,49, 92,128,148,141,111,64,0,-64,-111,-141,-148,-128,- 92,-49 ];
    y_coord= [ ...
      123,116, 91,42,-4,-67,-105,-119,-108,-119,-105,-67,-4,42, 91,116;
      129,132,112,62, 3,-57,-101,-110,-107,-110,-101,-57, 3,62,112,132;
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

  NODE(1,:)= NODE(1,:)*sqrt(2)/sqrt(1+ellip^2);
  NODE(2,:)= NODE(2,:)*sqrt(2)*ellip/sqrt(1+ellip^2);

  newmdl= mdl;
  newmdl.nodes= NODE';
  newmdl.elems= ELEM';
