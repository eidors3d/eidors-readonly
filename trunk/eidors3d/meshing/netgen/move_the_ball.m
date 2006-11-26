% Takes a netgen geo file with a ball in and moves it around
% First make the mesh
%Inputs are: (and defaults) 
%    electrodes_per_plane = 16;
%    number_of_planes     = 2;
%    movement_pattern     = 'spiral',
%                           'spirograph'
%                           'radial_move'
%                           'vertical_move'
%
%  Example
%    electrodes_per_plane = 16;  number_of_planes = 2; move_the_ball
%  Output:
%    output is saved to netgen_moving_ball.mat
%    
% (C) 2005 Bill Lionheart. Licensed under GPL v2
% - mods by Andy Adler to allow higher density electrode models
% $Id: move_the_ball.m,v 1.22 2006-11-26 00:59:56 aadler Exp $

% user input
if ~exist('electrodes_per_plane');  electrodes_per_plane= 16;    end
if ~exist('number_of_planes');      number_of_planes= 2;         end
if ~exist('movement_pattern');      movement_pattern= 'spiral';  end

% control mesh refinement: options are '-veryfine'; '-fine'; '';
if ~exist('finelevel');             finelevel= '';               end

fname =['tank_with_tank_',movement_pattern];


if ~exist('refine_electrodes');     refine_electrodes= 10;        end
if ~exist('tank_radius');           tank_radius= 15;              end
if ~exist('tank_height');           tank_height= 30;              end
if ~exist('electrode_width');       electrode_width = 0.5;        end
if ~exist('electrode_height');      electrode_height= 0.5;        end
if ~exist('rect_or_circ_electrode');rect_or_circ_electrode= 'C';  end

% number of simulations to do;
fno_max= 100;

first_plane_starts= tank_height/(number_of_planes+1);
height_between_centres = first_plane_starts;

[jnk,centres] = create_tank_mesh_ng( ...
   tank_radius, tank_height, ...
   rect_or_circ_electrode, ...
   log2(electrodes_per_plane), ...
   number_of_planes, ...
   first_plane_starts, ...
   height_between_centres, ...
   electrode_width, electrode_height, ...
   [fname,'0000'], ...
   refine_electrodes  );

% The msz file created here can be reused later
msz_file= [fname,'0000.msz'];

if ~isempty(finelevel);
   call_netgen([fname,'0000.geo'],[fname,'0000.vol'],msz_file, finelevel);
end


% Create Homog
stim_pat= mk_stim_patterns(electrodes_per_plane, number_of_planes, ...
              '{ad}','{ad}',{'meas_current'});
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
fmdl_save.elems=uint32(fmdl.elems);
      


% Load file and find places to modify
fid=fopen([fname,'0000.geo'],'r');
geo_homg= setstr(fread(fid)');
fclose( fid);
posn= findstr(geo_homg,'algebraic3d'); posn1=posn(1);
posn= findstr(geo_homg,'and  cyl;');   posn2=posn(1);

for fno= 1:fno_max
   % ensure memory isn't completely full
   eidors_cache clear all;

   r = 1.5;
   f_frac= fno/fno_max;
   switch movement_pattern
      case 'spiral'
         t=2*pi*f_frac * 4;
         x= tank_radius*f_frac*.7.*sin(t);
         y= tank_radius*f_frac*.7.*cos(t);
         z= f_frac*.7*tank_height + .15;

      case 'spirograph'
         % Parametric eqn from
         % online.redwoods.cc.ca.us/instruct/darnold/CalcProj/Fall98/CraigA/project3.htm
         t=2*pi*f_frac * 5;
%        t=2*pi*f_frac ; %  SLOW
         a=tank_radius;b=a*5/8;h=b*.8;
         x=(a-b)*cos(t) + h*cos( (a-b)/b*t );
         y=(a-b)*sin(t) - h*sin( (a-b)/b*t );
         z=tank_height/2;

      case 'vertical_move'
         t=  0;
         x= tank_radius*.5.*sin(t);
         y= tank_radius*.5.*cos(t);
         z= f_frac*(tank_height-4*r) + 2*r;

      case 'radial_move'
         t= 0;
         Rad_dist= (tank_radius-2*r)*f_frac;
         x= Rad_dist*sin(t);
         y= Rad_dist*cos(t);
         z= tank_height/2;


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
   call_netgen([fname_,'.geo'],[fname_,'.vol'],msz_file, finelevel);

   [fmdl,mat_idxs]= ng_mk_fwd_model( ...
     [fname_,'.vol'], centres, [], stim_pat);
   ball= mat_idxs{2};

   if 0 % Create a movie of the moving ball
      show_fem(fmdl);
      crop_model([],  inline('z>=30','x','y','z'))
      view(0,60);
      print('-dpng','-r100',[fname_,'.png']);
   end

   n_elem= size(fmdl.elems,1);
   elem_data= ones(n_elem,1);
   img=eidors_obj('image','netgen_problem', ...
                  'fwd_model',fmdl, 'elem_data', elem_data);
   img.elem_data(ball)= 2;
   vi(fno)= fwd_solve( img);
   vi(fno).time = fno;
end

save_filename= ...
   sprintf('tr=%f_th=%d_E=%s_EP=%d_NP_%d_EW=%f_EH=%f_RE=%d', ...
   tank_radius, tank_height, ...
   rect_or_circ_electrode, ...
   log2(electrodes_per_plane), ...
   number_of_planes, ...
   first_plane_starts, ...
   height_between_centres, ...
   electrode_width, electrode_height, ...
   [fname,'0000'], ...
   refine_electrodes  );
save(save_filename,'vi','vh','fmdl_save');


% Example to reconstruct images 
%imdl= mk_common_model('n3r2');
%imdl.fwd_model.stimulation= stim_pat;
%img= inv_solve(imdl,vi,vh);


