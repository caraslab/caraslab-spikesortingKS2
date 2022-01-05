function caraslab_kilosort(Savedir, rootH)
%This function runs kilosort2.
%
%Input variables:
%
%       Savedir: path to folder containing data directories. Each directory
%                should contain a binary (-dat) data file and
%                a kilosort configuration (config.mat) file. 
%
%       rootH:  path to Kilosort temp file; better if a fast SSD

%Written by ML Caras March 26, 2019 
% Patched by M Macedo-Lima 9/8/20

% Modified from master_kilosort.m in kilosort2 master repo.


%Prompt user to select folder
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end


%For each data folder...
for i = 1:numel(datafolders)
    clear ops rez
    close all

    cur_path.name = datafolders{i};
    cur_savedir = fullfile(Savedir, cur_path.name);

    %Load in configuration file (contains ops struct)
    % Catch error if -mat file is not found and skips folder
    try
        load(fullfile(cur_savedir, 'config.mat'));
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            fprintf('\n-mat file not found\n')
            continue
        else
            fprintf(ME.identifier)
            fprintf(ME.message)
            break
        end
    end

    % Sometimes kilosort will fail because it doesn't detect enough
%         % spikes in a batch. The batch size is determined by ops.NT
    % Increase ops.NT and rerun kilosort
    max_kilosort_attempts = 5; % Try this iteratively 5 times
    cur_attempts = 1;
    while cur_attempts <= max_kilosort_attempts
        try
            % Check if binary file exists in specified path, otherwise try to find it
            % in current folder and change paths appropriately
            if ~isfile(ops.fbinary)
                ops.fbinary = fullfile(cur_savedir, dir(fullfile(cur_savedir, '*CLEAN.dat')).name);
                if ~isfile(ops.fbinary)
                    fprintf(['Could not find binary file in: \n' cur_savedir]);
                    break
                else
                    ops.fclean = ops.fbinary;
                    save(fullfile(cur_savedir, 'config.mat'), 'ops');
                end
            end

            call_kilosort(cur_savedir, ops, rootH)
            % kilosort successful. Move on to next file

        catch ME
            % catch error resulting from few spikes per batch and rerun with
            % bigger ops.NT
            fprintf(['Kilosort failed with file:\n' cur_savedir]);
            fprintf('\nThe following error was produced: ')
            fprintf([ME.identifier '\n'])
            fprintf([ME.message '\n'])

            % these are the exceptions that can be generated when too few spikes
            % are found in a batch. Not sure why so many different
            % ones...
            if strcmp(ME.identifier, 'MATLAB:UndefinedFunction') || ...
                    strcmp(ME.identifier, 'MATLAB:eig:matrixWithNaNInf') || ...
                    strcmp(ME.identifier, 'parallel:gpu:kernel:LaunchFailure')

                cur_attempts = cur_attempts + 1;
                if cur_attempts > max_kilosort_attempts
                    % Failed 5 times, skipping file
                    fprintf(['\nKilosort failed 5x with file:\n' cur_savedir])
                    fprintf(['\nSkipping..........'])
                    break
                end

                % Double NT size then re-run
                ops.NT = ops.NT*2;
                save(fullfile(cur_savedir, 'config.mat'), 'ops', '-v7.3');

                fprintf('\nIncreasing batch size and trying again... Attempt: %d\n', cur_attempts)

                % Reload ops
                load(fullfile(cur_savedir, 'config.mat'));
                continue
            else
                fprintf('The following error was produced: ')
                fprintf(ME.identifier)  % print unknown exception
                rethrow(ME)
            end
        end
        break
    end

end


function call_kilosort(cur_savedir, ops, rootH)
    %Start timer
    t0 = tic;
    rootZ = cur_savedir; % the raw data binary file is in this folder

%         ops.trange = [180 Inf]; % time range to sort

    ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD

    % MML edit: change fbinary to be the clean, filtered file. I also
    % put a flag in the Kilosort gpufilter.m function to decide whether
    % to filter during Kilosort. Default is ops.kilosort_filter=0
    % because I filter it before.
    ops.fbinary = ops.fclean;

    %% this block runs all the steps of the algorithm
    fprintf('Looking for data inside %s \n', rootZ)

    rez = preprocessDataSub(ops);

    % time-reordering as a function of drift
    rez = clusterSingleBatches(rez);

    % saving here is a good idea, because the rest can be resumed after loading rez
    save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

    % main tracking and template matching algorithm
    rez = learnAndSolve8b(rez);

    % final merges
    rez = find_merges(rez, 1);

    % final splits by SVD
    rez = splitAllClusters(rez, 1);

    % final splits by amplitudes
    rez = splitAllClusters(rez, 0);

    % decide on cutoff
    rez = set_cutoff(rez);

    fprintf('found %d good units \n', sum(rez.good>0))

    % write to Phy
    fprintf('Saving results to Phy  \n')

    rezToPhy(rez, rootZ);

    tEnd= toc(t0);
    fprintf('Kilosort done in: %d minutes and %f seconds\n', floor(tEnd/60),rem(tEnd,60));

    delete(ops.fproc);






