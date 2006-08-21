% Takes a netgen geo file with a ball in and moves it around
% First make the mesh
% (C) 2005 Bill Lionheart. Licensed under GPL v2
% - mods by Andy Adler to allow higher density electrode models
% $Id: move_the_ball.m,v 1.15 2006-08-21 23:19:14 aadler Exp $

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
