function [d_insp,d_expi,W] = get_vent_meas( nn, imdl, ROI, FRATE, fidb);
   if nn(1:3)=='SS1'  % simulations
      idx= str2num( nn(11) );
      if idx==1;
         load montreal_data_1995;
         dd= [zc_resp, zc_resp, zc_resp, zc_resp, zc_resp];
         W = calc_reciproc_error(imdl, double( zc_resp ));
         FRATE=7;
      else
         load sim_radmove
         d_expi= vh;
         d_insp= -100*(vi(:, 105+idx)-vh)+vh;
         W = speye(length(vh));
         return;
      end
   else
      [dir, fname, ext] = fileparts(nn);
      ename = fullfile( '../breaths/' ,...
                        sprintf('%s-expi.mat',fname)        );
      iname = fullfile( '../breaths/' ,...
                        sprintf('%s-insp.mat',fname)        );
      if exist(ename, 'file') && exist(iname, 'file')
         fprintf('Using pre-calculated events for %s\n', fname);
         jnk = load(ename, 'dsv');
         d_expi = jnk.dsv;
         jnk = load(iname, 'dsv');
         d_insp = jnk.dsv;
         W = speye(length(d_insp));
         return
      end
      
      dd= eidors_readdata(nn);
      W = calc_reciproc_error(imdl, dd);
   end

   ok=1; % no problems with FRC indentification
   window= -0:0; owindow = ones(1,length(window));

   remove_pts = [];
   switch nn(end+(-14:0)) % special handling of files
      case 'S05/S05-004.get'
         dd= dd(:,1:725);
      case 'S05/S05-006.get'
         dd= dd(:,1:525);
      case 'S12/S12-006.get'
         remove_pts= linspace(550,605,4);
      case 'S16/S16-004.get'
         dd= dd(:,1:650);
   end
   
%     imdl.meas_icov = calc_reciproc_error(imdl,dd);
   ii = inv_solve( imdl, mean(dd,2), dd);   
   [einsp,eexpi] = find_frc( ii, ROI, FRATE, nn, ok, remove_pts);   
   fprintf('[%s: Look]',nn); pause

   f_expi= eexpi(:)*owindow + ones(length(eexpi),1)*window ;
   f_insp= einsp(:)*owindow + ones(length(einsp),1)*window ;
   d_expi= mean( dd(:,f_expi), 2);
   d_insp= mean( dd(:,f_insp), 2);
   
   
   if exist('ename','var')
      dsv = d_expi;
      save(ename,'dsv');
      dsv = d_insp;
      save(iname,'dsv');
   end
   
   
