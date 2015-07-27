function data = fwd_solve_peits(img, options)
%%fwd_solve_peits performs a forward solve in parallel (using PEITS)
%
% TODO: cache models and only update conductivity
% also do avoid repartitioning to save time!
%
% See also: https://users.dune-project.org/projects/dune-peits
%
%
% Usage Example:
%
% load('D:\Various\EITtoolbox\simulEITor\src\models\fbn\thoraxElecs_Coarse.mat')
% nSkip = 4
% [mdlThoraxElecs.stimulation, mdlThoraxElecs.meas_select] = mk_stim_patterns(32,1,[0,nSkip+1],[0,nSkip+1],{'no_rotate_meas', 'meas_current'},1);
% img = mk_image(mdlThoraxElecs,1);
% img.fwd_model.electrode = img.fwd_model.electrode(1:32)
% vtg = fwd_solve_peits(img);
%
%%
% Author: Fabian Braun <fbn@csem.ch>, September 2014
%%

%% unit testing?
if (nargin == 1) && (strcmpi(img, 'unit_test'))
    run_unit_test();
    return;
end

%% parse and default the options
if ~exist('options', 'var')
  options = [];
end
options = parse_and_default_options(options);



%% path settings
% where is PEITS located?
PEITS_PATH = getenv('PEITS_PATH');
if isempty(PEITS_PATH)
  % TODO: make this a bit more sophisticated
  % check for existance of peits and allow for reentering the path if not correct, etc.
  PEITS_PATH = input('Please enter the path to PEITS: ', 's');
  setenv('PEITS_PATH', PEITS_PATH);  
end
warning('move PEITS-specific m-files to PEITS directory!')
% addpath([PEITS_PATH, filesep, 'matlab']);

% the data dir where we store our files
FilePath = [PEITS_PATH, filesep, 'data', filesep];
if ~exist(FilePath, 'dir')
  mkdir(FilePath);
end


%% do caching in order to avoid regeneration of all files
CacheDependency = { img.fwd_model };
CacheName = ['fwd_solve_peits'];
[CacheOut, CacheId] = eidors_obj('get-cache', CacheDependency, CacheName);

if isempty(CacheOut)
    % set a dummy cache value to get an ID, 
    % actually here the eidors cache is just used to get a unique ID
    CacheId = eidors_obj('set-cache', CacheDependency, CacheName, 'foobar');
end

% generate various file names to either a) generate or b) reload files
CurrProtFileName = [FilePath, filesep, 'current_protocol_', CacheId, '.txt'];
dgfFileName = ['mesh_', CacheId, '.dgf'];
ElecFileName = [FilePath, filesep, 'electrode_nodes_', dgfFileName(1:end-4), '.txt'];
dgfFilePath = FilePath;
SigmaFileName = [FilePath, filesep, 'sigma_', CacheId, '.bin'];

if isempty(CacheOut)
    %% a pity: we'll have to generate all files since nothing was cached before

    %% somewhen do delete all unused files but keep ones in cache
    if options.cleanup_files
      cleanup_peits_files(PEITS_PATH);    
      % this is actually not such a good idea as it will delete all other cached files,
      % we should rather check which ones are not in cache anymore and delete only those!
    end

    %% parse stimulation patterns and create protocol txt file
    nElectrodes = length(img.fwd_model.stimulation);    
    InjCurrent = nan(1,nElectrodes);
    
    for iElec = 1:nElectrodes
        % go through all the different current injection pairs
        InjPos = find(img.fwd_model.stimulation(iElec).stim_pattern > 0); % current source
        InjNeg = find(img.fwd_model.stimulation(iElec).stim_pattern < 0); % current sink
        assert((length(InjPos) == 1) && (length(InjNeg) == 1), ...
            'only a pair of injecting electrodes supported!');
        
        InjCurrent(iElec) = img.fwd_model.stimulation(iElec).stim_pattern(InjPos);
        assert(InjCurrent(iElec) == -img.fwd_model.stimulation(iElec).stim_pattern(InjNeg));
        
        for iMeas = 1:size(img.fwd_model.stimulation(iElec).meas_pattern,1)
            % go through all the different voltage measurement pairs
            MeasPos = find(img.fwd_model.stimulation(iElec).meas_pattern(iMeas,:) == 1);
            MeasNeg = find(img.fwd_model.stimulation(iElec).meas_pattern(iMeas,:) == -1);
            assert((length(MeasPos) == 1) && (length(MeasNeg) == 1), ...
                'only a pair of injecting electrodes supported!');
            
            ActualPattern = [InjPos,InjNeg,MeasPos,MeasNeg];
            
            if (iMeas == 1) && (iElec == 1)
                dlmwrite(CurrProtFileName, ActualPattern);
            else
                dlmwrite(CurrProtFileName, ActualPattern, '-append');
            end
        end
    end
    
    assert(numel(unique([img.fwd_model.electrode.z_contact])) == 1, 'only one common contact impedance supported');
    assert(numel(unique(InjCurrent)) == 1, 'only one common injection current supported');    
    
    %% parse electrode nodes and generate node file   
    assert(nElectrodes == length(img.fwd_model.electrode), ...
        'Number of electrodes in stimulation and model do not agree!');
    
    for iElec=1:nElectrodes
        ElecNodes = sort(img.fwd_model.electrode(iElec).nodes);
        if iElec==1
            dlmwrite(ElecFileName, ElecNodes);
        else
            dlmwrite(ElecFileName, ElecNodes, '-append');
        end
    end
    
    %% export dgf file    
    dune_exporter(img.fwd_model.nodes, img.fwd_model.elems, img.elem_data, ...
        dgfFilePath, dgfFileName, ...
        zeros(nElectrodes, 3), ...
        img.fwd_model.nodes(img.fwd_model.gnd_node,:))
        
    %% conductivity is stored in dgf file directly, don't load separately
    SigmaFileName = [];    
else
    %% YES! we can rely on cached variables
    
    %% create a separate file for conductivity values
    SigmaIds = CacheOut;    
    save_sigma_vector_binary(SigmaFileName, img.elem_data(SigmaIds), SigmaIds);
    
    InjCurrent = unique(cell2mat(cellfun(@(x) abs(full(x)), ...
                    {img.fwd_model.stimulation.stim_pattern}, 'uniformoutput', false)));
    InjCurrent(InjCurrent == 0) = [];
                
    assert(length(InjCurrent) == 1, 'only one common injection current supported');
end


%% run the thingy here
% prepare settings
forward_settings = set_forward_default_values();
forward_settings.path = PEITS_PATH;
forward_settings.mesh.name = dgfFileName(1:end-4);
[~, CurrProtFileNameOnly, CurrProtFileExtOnly] = fileparts(CurrProtFileName);
forward_settings.protocol = [CurrProtFileNameOnly, CurrProtFileExtOnly];
forward_settings.contact_impedance = unique([img.fwd_model.electrode.z_contact]);
forward_settings.current = unique(InjCurrent);
forward_settings.do_jacobian = 0;
forward_settings.measORall = 1;
forward_settings.do_elec_volts = 1;

if ~isempty(CacheOut)    
    % in case of using the cached model we have to provide the updated conductivity (sigma) file 
    forward_settings.fem.io.load_sigma_separately = 1;
    [~, SigmaFileNameOnly, SigmaFileExtOnly] = fileparts(SigmaFileName);
    forward_settings.fem.io.separate_sigma_file = ['..', filesep, 'data', filesep, SigmaFileNameOnly, SigmaFileExtOnly];     
    display('using cached model with separate conductivities...');
else
    forward_settings.fem.io.load_sigma_separately = 0;    
    display('no cached model...');
end

%% finally, run the solver
[vtg,j,output,SigmaIds] = run_forward_solver(forward_settings, options.n_processes);

if isempty(CacheOut)
    % set a dummy cache value to get an ID, 
    % actually here the eidors cache is just used to get a unique ID
    CacheIdNew = eidors_obj('set-cache', CacheDependency, CacheName, SigmaIds);    
    assert(strcmp(CacheId,CacheIdNew));
end


%% create a data structure to return
% TODO: make this proper similar to what is done in: fwd_solve_1st_order
data.meas = vtg;
% data.meas= meas_from_v_els(v_els, fwd_model.stimulation);
data.time = NaN; % unknown
data.name = 'solved by fwd_solve_peits';
data.type = 'data';
% try; if img.fwd_solve.get_all_meas == 1
%    data.volt = v(1:pp.n_node,:); % but not on CEM nodes
% end; end
% try; if img.fwd_solve.get_all_nodes== 1
%    data.volt = v;                % all, including CEM nodes
% end; end

end

function options = parse_and_default_options(options)

  if isempty(options)
    options = struct;
  end

  if ~isfield(options, 'n_processes')
    if isunix()
        options.n_processes = 16; % UNIX	
    else
        options.n_processes = 4;  % Windows
    end
  end

  if ~isfield(options, 'cleanup_files')
    % before creating files for a new model, shall we delete all old files?
    options.cleanup_files = false;
  end
  
  % TODO: add more options here?!?
    
end

function run_unit_test


%% compare EIDORS and PEITS on a simple model
% eidors_debug('on', 'call_netgen')
fmdl = ng_mk_cyl_models(2,[32,1],[0.05]); 
nSkip = 4
[fmdl.stimulation, fmdl.meas_select] = mk_stim_patterns(32,1,[0,nSkip+1],[0,nSkip+1],{'no_rotate_meas', 'meas_current'},1);
img = mk_image(fmdl,1);
img.elem_data = rand(size(img.elem_data));
img.fwd_model.electrode = img.fwd_model.electrode(1:32);
meas = fwd_solve_peits(img);
meas2 = fwd_solve(img);

figure, plot(meas2.meas)
% hold on, plot(meas.meas*0.5659/1.235, 'r')  % why should they differ by such a factor?
hold on, plot(meas.meas, 'r')
legend('EIDORS', 'PEITS');
display('visually inspect the two U-shapes: stopping here and waiting for user');
keyboard;

% TODO: compare the two results programatically

% TODO: implement lots of other tests

%% TODO: add much more here
warning('improve unit_test!');

%% all done?
eidors_msg('@@@ unit test completed',1);

end

function cleanup_peits_files(PEITS_PATH)

	delete([PEITS_PATH, filesep, 'partitions', filesep, 'eidors*']);
	delete([PEITS_PATH, filesep, 'output', filesep, '*.bin']);
	delete([PEITS_PATH, filesep, 'data', filesep, '*eidors*']);
	delete([PEITS_PATH, filesep, 'data', filesep, 'current_protocol_id_*']);
	delete([PEITS_PATH, filesep, 'data', filesep, '*mesh_id_*']);
	delete([PEITS_PATH, filesep, 'data', filesep, 'electrode_nodes_id_*']);
	delete([PEITS_PATH, filesep, 'data', filesep, 'electrode_positions*_id_*']);
	delete([PEITS_PATH, filesep, 'data', filesep, 'sigma_id_*']);

end
