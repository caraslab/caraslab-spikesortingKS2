function compile_data_for_analyses(Savedir)
%
% This function loops through recording folders and extracts all relevant
% files into a folder called Data inside the parent directory. The purpose
% is to centralize all subjects' data into common directories
%


%Prompt user to select folder
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end


% Create data analysis folders in parent directory
full_path_split = split(Savedir, "/");
parent_path = strjoin(full_path_split(1:end-1), "/");
mkdir(fullfile(parent_path, 'Data'));
mkdir(fullfile(parent_path, 'Data', 'Behavioral performance'));
mkdir(fullfile(parent_path, 'Data', 'Breakpoints'));
mkdir(fullfile(parent_path, 'Data', 'Key files'));
mkdir(fullfile(parent_path, 'Data', 'Output'));
mkdir(fullfile(parent_path, 'Data', 'Quality metrics'));
mkdir(fullfile(parent_path, 'Data', 'ShankWaveform plots'));
mkdir(fullfile(parent_path, 'Data', 'Spike times'));
mkdir(fullfile(parent_path, 'Data', 'Unit database'));
mkdir(fullfile(parent_path, 'Data', 'Waveform measurements'));
mkdir(fullfile(parent_path, 'Data', 'DAC files'));

%For each data folder...
for i = 1:numel(datafolders)
        cur_path.name = datafolders{i};
        cur_savedir = fullfile(Savedir, cur_path.name);

        % Copy breakpoints
        filedirs = dir(fullfile(cur_savedir, 'CSV files', '*breakpoints.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'Breakpoints', cur_filedir.name));
            end
        end
        
        % Copy Key files
        filedirs = dir(fullfile(cur_savedir, 'CSV files', '*trialInfo.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'Key files', cur_filedir.name));
            end
        end
        filedirs = dir(fullfile(cur_savedir, 'CSV files', '*spoutTimestamps.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'Key files', cur_filedir.name));
            end
        end
        
        % Copy quality metrics
        filedirs = dir(fullfile(cur_savedir, 'CSV files', '*quality_metrics.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'Quality metrics', cur_filedir.name));
            end
        end
        
        % Copy ShankWaveform plots
        filedirs = dir(fullfile(cur_savedir, '*shankWaveforms*.pdf'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'ShankWaveform plots', cur_filedir.name));
            end
        end        
        
        % Copy Spike times
        filedirs = dir(fullfile(cur_savedir, '*cluster*.txt'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'Spike times', cur_filedir.name));
            end
        end      
        
        % Copy waveform measurements
        filedirs = dir(fullfile(cur_savedir, 'CSV files', '*waveform_measurements.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'Waveform measurements', cur_filedir.name));
            end
        end
        
        % Copy DAC files
        filedirs = dir(fullfile(cur_savedir, 'CSV files', '*DAC*.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(parent_path, 'Data', 'DAC files', cur_filedir.name));
            end
        end
end

% Lastly copy behavior files
filedirs = caraslab_lsdir(fullfile(Savedir, 'Behavior'));
filedirs = {filedirs.name};
for file_idx=1:length(filedirs)
    if ~isempty(filedirs)
        cur_filedir = fullfile(Savedir, 'Behavior', filedirs{file_idx});
        copyfile(cur_filedir, fullfile(parent_path, 'Data', 'Behavioral performance', filedirs{file_idx}));
    end
end