% Takes a netgen geo file with a ball in and moves it around
% First make the mesh
% (C) 2005 Bill Lionheart. Licensed under GPL v2
% - mods by Andy Adler to allow higher density electrode models
% $Id: move_the_ball.m,v 1.14 2006-08-21 23:16:39 aadler Exp $

fname ='tank_for_kalman_test_ball_';

refine_electrodes= 50;
tank_radius= 15;
tank_height= 30;
[tank_mdl1,centres] = create_tank_mesh_ng( ...
                        tank_radius, ...
                        tank_height, ...
                        'R', ... % rect_electrodes
                        4,  ... % log2_electrodes_per_plane,
                        2, ...  % no_of_planes,
                        tank_height/3, ... % first_plane_starts
                        10, ... % height_between_centres
                        1.5, ... % electrode_width
                        1.5, ... % electrode_height
                        [fname,'0000'], ...
                        refine_electrodes  );


% Create Homog
stim_pat= mk_stim_patterns(16,2,'{ad}','{ad}',{'meas_current'});
[fmdl,mat_idxs]= ng_mk_fwd_model( ...
  [fname,'0000.vol'], centres, [], stim_pat);
n_elem= size(fmdl.elems,1);
elem_data= ones(n_elem,1);
img=eidors_obj('image','netgen_problem', ...
               'fwd_model',fmdl, 'elem_data', elem_data);
vh= fwd_solve( img);
vh.time= 0;

if 0 % TEST CODE - check electrodes
   show_fem( fmdl);

   % verify np_fwd_parameters identifies
   % the correct elements
   pp= np_fwd_parameters( fmdl);
   for i= 1:16
      ee= reshape( pp.elec(i,:), 3, []);
      ee(:, all(ee==0,1)) = [];

      Xs = reshape(fmdl.nodes(ee,1),3,[]);
      Ys = reshape(fmdl.nodes(ee,2),3,[]);
      Zs = reshape(fmdl.nodes(ee,3),3,[]);
      h=patch(Xs,Ys,Zs, 'b');
      set(h, 'FaceLighting','none', 'CDataMapping', 'direct' );

      elec_area= 0;
      for j= 1: size(ee,2)
         elec_area = elec_area + triarea3d( fmdl.nodes( ee(:,j),:) );
      end

      fprintf('elec#%d: area=%f\n',i,elec_area);pause;
   end
end

fmdl_save=fmdl;
% save memory
fmdl_save.nodes=uint16(fmdl.nodes);
fmdl_save.elems=uint16(fmdl.elems);
      


% Load file and find places to modify
fid=fopen([fname,'0000.geo'],'r');
geo_homg= setstr(fread(fid)');
fclose( fid);
posn= findstr(geo_homg,'algebraic3d'); posn1=posn(1);
posn= findstr(geo_homg,'and  cyl;');   posn2=posn(1);

fno_max= 100;
% go around circle and move to boundary
for fno= 1:fno_max
   % ensure memory isn't completely full
   eidors_cache clear all;

   r = 1.5;
   f_frac= fno/fno_max;
   switch 'spiral'
      case 'spiral'
         t=2*pi*f_frac * 4;
         x= tank_radius*f_frac*.7.*sin(t);
         y= tank_radius*f_frac*.7.*cos(t);
         z= f_frac*.7*tank_height + .15;

      case 'spirograph'
         % Parametric eqn from
         % online.redwoods.cc.ca.us/instruct/darnold/CalcProj/Fall98/CraigA/project3.htm
         t=2*pi*f_frac * 5; % GCD of a and b
         a=tank_radius;b=a*5/8;h=b*.8;
         x=(a-b)*cos(t) + h*cos( (a-b)/b*t );
         y=(a-b)*sin(t) - h*sin( (a-b)/b*t );
         z=tank_height/2;

      otherwise
         error('don''t understand movement options')
   end
   % Matlab's stupid inconsistency with \n is _really_ frunstrating!!
   geo_file= [ ...
       geo_homg(1:(posn1+11)), ...
       sprintf('solid ball = sphere(%3.3f,%3.3f,%3.3f;%3.3f);', ...
          x,y,z,r), ...
       10,geo_homg((posn1+12):(posn2+7)), ...
       ' and not ball',geo_homg((posn2+8):end), ...
       10,'tlo ball   -col=[0,1,0];', 10];
   fname_ = sprintf('%s%04d', fname, fno );
   fid=fopen([fname_,'.geo'],'w');
   fwrite(fid,geo_file);
   fclose( fid);
   call_netgen([fname_,'.geo'],[fname_,'.vol']);

   [fmdl,mat_idxs]= ng_mk_fwd_model( ...
     [fname_,'.vol'], centres, [], stim_pat);
   ball= mat_idxs{2};

   if 0 % Create a movie of the moving ball
      crop_model([],  inline('z>=30','x','y','z'))
      view(0,60);
      print('-dpng','-r100',[fname_,'.png']);
   end

   n_elem= size(fmdl.elems,1);
   elem_data= ones(n_elem,1);
   img=eidors_obj('image','netgen_problem', ...
                  'fwd_model',fmdl, 'elem_data', elem_data);
   img.elem_data(ball)= 2;
   vi{fno}= fwd_solve( img);
   vi{fno}.time = fno;
end

save netgen_moving_ball vi vh fmdl_save

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnstem ='tank_for_kalman_test_ball_';
fnstemnb=[fnstem,'no_ball'];
fnnb =[fnstemnb,'.geo'];

fn00= [fnstem,'00.geo'];

tank_radius =15;
tank_height=30;
CorR = 'C';  %C for circular R for rectangular
log2_electrodes_per_plane =4;
no_of_planes=2;
first_plane_starts=10; 
height_between_centres=10;
%electrode_width=1.5;
%electrode_height=2.5;
electrode_width=1.5;
electrode_height=1.5;

elecs_per_plane=2^log2_electrodes_per_plane;
nelec=no_of_planes*elecs_per_plane;


%Make list of centres
kel = 0;
for l=1:no_of_planes
 z = first_plane_starts + (l-1)*height_between_centres;
 for th =0:2*pi/elecs_per_plane:2*pi*(elecs_per_plane-1)/elecs_per_plane
    kel=kel+1;
    [x,y]=pol2cart(th+0.02*randn(1),tank_radius);  %jiggle to avoid netgen errors
    dirn = [x,y,0];
    centres(kel,:)= [x,y,z]; % keep the centres
    
 end;
end;


tank_mdl1 = create_tank_mesh_ng_centres( ...
                        tank_radius, tank_height, ...
                        CorR,centres,...
                        electrode_width,electrode_height,fnstemnb);
% comment out if already made
r = tank_radius/10; %that is radius of ball
%Now put a ball in somewhere
x=0;y=0;z=tank_height/2;
%fnnb
fidnb=fopen(fnnb);
fidout =fopen(fn00,'w');
foundit=0;
  while 1
    tline = fgetl(fidnb);
    if ~ischar(tline), break, end
    if ~isempty( findstr(tline,'algebraic3d')) & ~foundit
         fprintf(fidout,'%s\n',tline);
         foundit = 1;
         tline=   sprintf('solid ball = sphere(%3.3f,%3.3f,%3.3f;%3.3f);',x,y,z,r);
    end
%    disp(tline)
% take the ball from the tank
if strfind(tline,'and  cyl;')  %care  - two spaces!
    disp(tline);
    tline='and cyl and not ball;'
end
fprintf(fidout,'%s\n',tline);
  end  %while
fprintf(fidout,' tlo ball   -col=[0,1,0];\n')
fclose(fidout);
fclose(fidnb);

%Made the base file
fnin = fn00;

ntimes = 4;% 32;

for itime = 1:ntimes
  fnouts=sprintf('tank_for_kalman_test_ball_%02d',itime);
  fnout=[fnouts,'.geo'];
  fidout=fopen(fnout,'w');

  fidin=fopen(fn00,'r');
  fidout=fopen(fnout,'w');

  % Ball starts at center on bottom and
  % Rotates towards the top and outward
  z = (itime/(ntimes+1))*(tank_height - 2*r)+ r;
  R = ((itime-1)/ntimes)*(tank_radius-r);
  th = itime*2*pi/(ntimes/4);
  [x,y]=pol2cart(th,R);
  % Round to 3dp as netgen hates -0.000
  x= round(x*1000)/1000;
  y= round(y*1000)/1000;
  z= round(z*1000)/1000;
  r= round(r*1000)/1000;

  while 1
    tline = fgetl(fidin);
    if ~ischar(tline), break, end
    if ~isempty( findstr(tline,'solid'))
          if  ~isempty( findstr(tline,'ball'))
              tline=   sprintf('solid ball = sphere(%3.3f,%3.3f,%3.3f;%3.3f);',x,y,z,r);
          end
    end
%    disp(tline)
  if strfind(tline,  'and  cyl;')
      tline = 'and  cyl and not ball;'
  end    
    fprintf(fidout,'%s\n',tline);
  end  %while
  fclose(fidout);
  fclose(fidin);
  meshfn=[fnouts,'.vol'];

  call_netgen(fnout,meshfn);
end %for

% Now read them in again and do the calculations

for itime = 1:ntimes
  fnins=sprintf('tank_for_kalman_test_ball_%02d',itime);
  fnin=[fnins,'.vol'];
  [srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(fnin);
  the_mat_indices = unique(mat_ind);
  if size(the_mat_indices,1)~=2
     error('This mesh does not have two material indices');
  end
  %Assume the sphere is smaller than the rest of tank
  [junk,k]=min([length(find (mat_ind==the_mat_indices(1))), ...
                length(find (mat_ind==the_mat_indices(2)))]);
  mat_index_sphere= the_mat_indices(k);
  elements_in_sphere = find(mat_ind==mat_index_sphere);

  % Construct mesh object
  tank_mdls(itime).type= 'fwd_model';
  tank_mdls(itime).name= sprintf('Tank model with ball time %d',itime);
  tank_mdls(itime).nodes= vtx;
  tank_mdls(itime).elems= simp;
  tank_mdls(itime).boundary= srf;
  [elec,sels] = ng_tank_find_elec(srf,vtx,bc,centres);
  if size(elec,1) ~= nelec
     error('Failed to find all the electrodes')
  end
  for i=1:nelec
      electrodes(i).z_contact= 1.0;  %well you can change that later!
      electrodes(i).nodes=     unique( elec(i,:) );
  end
  perm_sym='{n}';
  tank_mdls(itime).type = 'fwd_model';
  tank_mdls(itime).gnd_node=           1;
  tank_mdls(itime).electrode =         electrodes;
  tank_mdls(itime).misc.perm_sym =     perm_sym;

  tank_mdls(itime).solve= 'np_fwd_solve';
  tank_mdls(itime).jacobian= 'np_calc_jacobian';
  tank_mdls(itime).system_mat= 'np_calc_system_mat';
  tank_mdls(itime).stimulation =mk_stim_patterns( ...
                     2^log2_electrodes_per_plane, ...
                     no_of_planes, '{op}', '{ad}',{'meas_current'});


  cond = 1.0*ones(size(simp,1),1);
  cond(elements_in_sphere) =  2.0;  % The contrast
  tank_img(itime) = eidors_obj('image', sprintf('ball image t=%d',itime), ...
                       'elem_data', cond, ...
                       'fwd_model', tank_mdls(itime) );
  tank_data(itime)=fwd_solve(tank_img(itime));

  if itime==1
     cond = 1.0*ones(size(simp,1),1);
     tank_img_homg= tank_img(itime);
     tank_img_homg.elem_data = cond;
     tank_data_homg=fwd_solve(tank_img_homg);
  end


end %for


save('all_models_for_moving_ball', ...
              'tank_img','tank_data', ...
              'tank_img_homg','tank_data_homg');


return;

% Now, create 2D simulations with 16 electrodes in a plane
sim_pat2d= mk_stim_patterns(16,1,'{ad}','{ad}',{'no_meas_current'});
th= tank_img_homg;
th.fwd_model.stimulation= sim_pat2d;
for i=1:nelec
  th.fwd_model.electrode(i).z_contact= 0.01;
end
th.fwd_model.electrode(17:32)=[]; %remove top layer
tank_data_2d_h = fwd_solve(th);
for t=1:ntimes
    th= tank_img(t);
    th.fwd_model.stimulation= sim_pat2d;
    for i=1:nelec
      th.fwd_model.electrode(i).z_contact= 0.01;
    end
    th.fwd_model.electrode(17:32)=[]; %remove top layer
    tank_data_2d(t) = fwd_solve(th);
end


