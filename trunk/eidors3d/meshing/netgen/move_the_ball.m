% Takes a netgen geo file with a ball in and moves it around
% First make the mesh
% (C) 2005 Bill Lionheart. Licensed under GPL v2
% $Id: move_the_ball.m,v 1.9 2006-07-29 15:38:01 aadler Exp $


fnstem ='tank_for_kalman_test_ball_';
fnstemnb=[fnstem,'no_ball'];
fnnb =[fnstemnb,'.geo'];

fn00= [fnstem,'00.geo'];

tank_radius =15;
tank_height=30;
CorR = 'R';
log2_electrodes_per_plane =4;
no_of_planes=2;
first_plane_starts=10; 
height_between_centres=10;
electrode_width=1.5;
electrode_height=2.5;
nelec=no_of_planes*2^log2_electrodes_per_plane;

[tank_mdl1,centres] = create_tank_mesh_ng( ...
                        tank_radius, tank_height, ...
                        CorR,log2_electrodes_per_plane, ...
                        no_of_planes,first_plane_starts, ...
                        height_between_centres, ...
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
    fprintf(fidout,'%s\n',tline);
  end  %while
fprintf(fidout,' tlo ball   -col=[0,1,0];\n')
fclose(fidout);
fclose(fidnb);

%Made the base file
fnin = fn00;

ntimes = 32;

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
    fprintf(fidout,'%s\n',tline);
  end  %while
fclose(fidout);
fclose(fidin);
meshfn=[fnouts,'.vol'];
status= system(sprintf( ...
        'ng -batchmode -geofile=%s  -meshfile=%s ',fnout,meshfn));
   
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

