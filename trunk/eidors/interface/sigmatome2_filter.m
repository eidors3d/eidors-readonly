function [Filter, stim_pattern]= sigmatome2_filter(test);
% SIGMATOME2_FILTER:  Hardware filter and stim_patterns for Sigmatome II device
%
% Usage:
%  [Filter, stim_pat]= sigmatome2_filter;
%    Filter   is a 416x208 Hardware filter
%    stim_pat is a EIDORS stim_pattern structure
%
% Example:
%   imdl= mk_common_model('c2t2',16);
%   [Filter,stim_pat] = sigmatome2_filter;
%   imdl.fwd_model.stimulation = stim_pat;
%   imdl.fwd_model.jacobian = @jacobian_filtered;
%   imdl.fwd_model.jacobian_filtered.jacobian = @jacobian_adjoint;
%   imdl.fwd_model.jacobian_filtered.filter   = Filter;
%   imdl.meas_icov = speye( size(Filter,1) );
%   imdl.fwd_model = rmfield(imdl.fwd_model, 'meas_select');
%
%   img= mk_image(imdl);
%   vh= fwd_solve(img);  vh = Filter * vh.meas;
%   img.elem_data(50)=1.1;
%   vi= fwd_solve(img);  vi = Filter * vi.meas;
%
%   imdl.hyperparameter.value = 0.03;
%   imgr= inv_solve(imdl, vh, vi);
%   show_slices(imgr);

% (C) 2015 Andy Adler. License: GPL version 2 or version 3
%   Based on information from Robert Guardo and Herve Gagnon
% $Id$

if nargin>=1 && strcmp(test,'UNIT_TEST'); do_unit_test; return; end

Filter   = calc_sigmatome2_filter;
protocol = ComputeAllConfigs([0 1 2 3], 16) + 1;
stim_pattern = stim_meas_list( protocol );


% Generate the filter function of the sigmatome system
function filter = calc_sigmatome2_filter
   Filter3= [0,         0.002992185737274, 0.002885352918496,-0.001992209356618, ...
     0.000247868161165, 0.001320969614758,-0.002994506380400, 0.000217287704644, ...
     0.000820760030266,-0.002388560206626, 0.005467093174149,-0.000465034175608, ...
    -0.014049400948852, 0.011968526958771, 0.017269336581277,-0.033327865383464, ...
    -0.000687996077078, 0.065884070897957,-0.055509942269519,-0.084217637095161, ...
     0.414790173689917, 0.947607669916644, 0.723749489957398, 0.088544295521622, ...
    -0.118281839343952, 0.032044694797133, 0.043642225555493,-0.035476390241121, ...
    -0.007655680300629, 0.022891001199580,-0.005076656322158,-0.007794874431280, ...
     0.006724761261718,-0.000064773581282,-0.002313251865051, 0.001234620138549, ...
    -0.002848224928827,-0.000772726894395, 0.002559973165365,-0.001297698809235, ...
     0.000569993353282, 0.002828133181959, 0.000426156082041];
   filter = sparse(416,208);
   filter(1:length(Filter3),3) = Filter3';
   for i=1:208; 
      filter(:,i) = circshift(filter(:,3),+2*(i-3));
   end

function [Configs Mesptr Curptr Index] = ComputeAllConfigs(InitConfig, N);
%
%ComputeAllConfigs: Cette fonction calcule toutes les configurations
%                  [Source, Puit, Suiveur, Inverseur] a partir de la
%                  configuration ititiale et du nombre d'electrodes.
%
%SYNTAXE:  Configs = ComputeAllConfigs(InitConfig, NbreElectrodes);
%
%INPUT:    InitConfig(1 x 4):  Configuration initiale de mesure exprimee dans 
%                              le format suivant : [Source, Puit, Suiveur, 
%                              Inverseur]. Les electrodes sont numerotees de
%                              0 a (NbreElectrodes-1).
%
%          N:                  Nombre d'electrodes.
%
%OUTPUT:   Configs(m x 4):     Configurations d'electrodes : [Source, Puit, 
%                              Suiveur, Inverseur] pour les "m" mesures prises
%                              par le scanhead.
%
%          Mesptr(2 x N):      Chaque colonne correspond aux positions
%                              d'electrodes qui forment une des "N" paires 
%                              d'electrodes utilises pour effectuer les mesures 
%                              de tension.
%
%          Curptr(2 x N):      Chaque colonne correspond aux numeros des noeuds
%                              qui forment une des "N" paires d'electrodes 
%                              utilises pour injecter les courants.
%
%          Index(m x 2):       Pour les "m" mesures prises par le scanhead,
%                              doublet indiquant le numero de la paire 
%                              d'electrodes qui sert a applique le courant et 
%                              le numero de la paire qui sert a effectuer la 
%                              mesure.

% Copyright (c) 2009 Hervé Gagnon, Ecole Polytechnique de Montréal.

   % Validation du parametre InitConfig
   if (size(InitConfig,1) ~= 1 || size(InitConfig,2) ~= 4)
       error('Le parametre "InitConfig" doit etre de dimension (1 x 4)!');
   end

   % Validation du parametre NbreElectrodes
   if (length(N) ~= 1)
       error('Le parametre "NbreElectrodes" doit etre un scalaire!');
   end

   SrcePosInit = InitConfig(1);
   SinkPosInit = InitConfig(2);
   FollPosInit = InitConfig(3);
   InvtPosInit = InitConfig(4);

   GapCourant = abs(SrcePosInit - SinkPosInit);
   GapTension = abs(FollPosInit - InvtPosInit);

   if ((GapCourant == N/2) && (GapTension == N/2))
       m = N*(N - 2);
   elseif((GapCourant == GapTension) || (GapCourant + GapTension == N))
       m = N*(N - 3);
   else
       m = N*(N - 4);
   end
      
   SrcePos = SrcePosInit;
   SinkPos = SinkPosInit;
   FollPos = FollPosInit;
   InvtPos = InvtPosInit;

   Configs = zeros(m, 4);

   for i = 1:m
       Configs(i,:) = [SrcePos SinkPos FollPos InvtPos];
       
       SrcePos = IncrementVarModuloN(SrcePos, N);
       SinkPos = IncrementVarModuloN(SinkPos, N);
       FollPos = IncrementVarModuloN(FollPos, N);
       InvtPos = IncrementVarModuloN(InvtPos, N);
       
       if ((SrcePos == SrcePosInit) && (SinkPos == SinkPosInit))
          FollPos = IncrementVarModuloN(FollPos, N);
          InvtPos = IncrementVarModuloN(InvtPos, N);
          while (FollPos == SrcePos || FollPos == SinkPos || InvtPos == SrcePos || InvtPos == SinkPos)
              FollPos = IncrementVarModuloN(FollPos, N);
              InvtPos = IncrementVarModuloN(InvtPos, N);
          end
       end
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%
   Curptr = zeros(2,N);
   Mesptr = zeros(2,N);

   Curptr(:,1) = [SrcePosInit; SinkPosInit];
   Mesptr(:,1) = [FollPosInit; InvtPosInit];

   for i = 2:N
       Curptr(1,i) = IncrementVarModuloN(Curptr(1,i-1), N);
       Curptr(2,i) = IncrementVarModuloN(Curptr(2,i-1), N);
       Mesptr(1,i) = IncrementVarModuloN(Mesptr(1,i-1), N);
       Mesptr(2,i) = IncrementVarModuloN(Mesptr(2,i-1), N);
   end

   Index = zeros(m, 2);
   for i = 1:m
       Index(i,1) = find((Configs(i,1) == Curptr(1,:)) & (Configs(i,2) == Curptr(2,:)));
       Index(i,2) = find((Configs(i,3) == Mesptr(1,:)) & (Configs(i,4) == Mesptr(2,:)));
   end

function VarOut = IncrementVarModuloN(VarIn, N);
   VarOut = VarIn + 1;
   if (VarOut == N)
       VarOut = 0;
   end

function do_unit_test
   imdl= mk_common_model('c2t2',16);
   [Filter,stim_pat] = sigmatome2_filter;
   imdl.fwd_model.stimulation = stim_pat;
   imdl.fwd_model.jacobian = @jacobian_filtered;
   imdl.fwd_model.jacobian_filtered.jacobian = @jacobian_adjoint;
   imdl.fwd_model.jacobian_filtered.filter   = Filter;
   imdl.meas_icov = speye( size(Filter,1) );
   imdl.fwd_model = rmfield(imdl.fwd_model, 'meas_select');

   img= mk_image(imdl);
   vh= fwd_solve(img);  vh = Filter * vh.meas;
   img.elem_data(50)=1.1;
   vi= fwd_solve(img);  vi = Filter * vi.meas;

   imdl.hyperparameter.value = 0.03;
   imgr= inv_solve(imdl, vh, vi);
   show_slices(imgr);

   imgr = calc_slices(imgr);
   max_imgr = max(imgr(:));
   unit_test_cmp('Max', max_imgr, 6.683411990674756e-04, 1e-8);
   unit_test_cmp('Loc', find(imgr>0.99*max_imgr),  ...
      [1755; 1756; 1818; 1819; 1820; 1882; 1883; 1946; 1947]);
