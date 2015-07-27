function [v,j,output,id] = run_forward_solver(forward_settings, processes)

if forward_settings.path(end) ~= filesep
    forward_settings.path(end+1) = filesep;
end

% write the settings to the parameter files
if ~exist([forward_settings.path, 'data', filesep, forward_settings.mesh.name, '.dgf'], 'file')
    error('This mesh has not yet been written to DGF format. Use dune_exporter.m to do that.');
end
if ~exist([forward_settings.path, 'partitions', filesep, forward_settings.mesh.name, '.dgf.', num2str(processes), '_0'], 'file')
    disp('This mesh has not yet been partitioned to the defined number of processes. This will take a while.');
    forward_settings.mesh.do_partitions = true;
end
if isa(forward_settings.protocol, 'double')
    disp('Creating new protocol file "data', filesep, 'protocol.txt".');
    fid = fopen([forward_settings.path,'data', filesep, 'protocol.txt'],'w');
    fprintf(fid,'%f,%f,%f,%f\n',forward_settings.protocol');
    fclose(fid);
    forward_settings.protocol = 'protocol.txt';
end

write_parameter_files(forward_settings);

% run the solver
if isunix()
	fwdSolverCmd = ['cd ', forward_settings.path, filesep, 'src/; mpirun -np ', num2str(processes), ' ./dune_peits; exit']
else
  % TODO: set the cygwin path as option and not hard-coded
  path_cygwin = strrep(forward_settings.path, filesep, '/');
  path_cygwin = strrep(path_cygwin, ':/', '/');
  path_cygwin = ['/cygdrive/', path_cygwin];
	fwdSolverCmd = ['"C:\cygwin64\bin\bash" -li -c "cd ', path_cygwin, '/src/; mpirun -np ', num2str(processes), ' ./dune_peits; exit"']
end
[status,output] = system(fwdSolverCmd, '-echo');

% check that everything went well
if status
    disp(output);
    error('Something went wrong when running the solver: stopping here!');
end

% read out the results
if forward_settings.do_elec_volts
    % find out which is the newest file and load it
    d = dir([forward_settings.path,filesep,'output',filesep,'electrodevoltages*.bin']);
    [~,index] = max([d.datenum]);
    vFileName = [forward_settings.path,filesep, 'output',filesep, d(index).name];
    display(['loading voltages from: ' vFileName]);
    v = load_electrode_voltages_binary(vFileName);
else
    v = 0;
end


d = dir([forward_settings.path,'output', filesep, 'sigmavector*.bin']);
[~,index] = max([d.datenum]);
[id,~] = load_sigma_vector_binary([forward_settings.path,'output', filesep, d(index).name]);


if forward_settings.do_jacobian
    % find out which is the newest file and load it. also, sort the
    % jacobian by the element ID's
    d = dir([forward_settings.path,'output', filesep, 'sigmavector*.bin']);
    [~,index] = max([d.datenum]);
    [id,~] = load_sigma_vector_binary([forward_settings.path,'output',filesep, d(index).name]);
    d = dir([forward_settings.path,'output', filesep, 'jacobian*.bin']);
    [~,index] = max([d.datenum]);
    j_unsorted = load_jacobian_binary([forward_settings.path,'output', filesep, d(index).name]);
    j = zeros(size(j_unsorted));
    j(:,id) = j_unsorted;
else
    j = 0;
end
