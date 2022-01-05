function [wf, gwfparams] = try_load_previous_gwfparams(day_dir, bp_filter, only_good)    
    % Load waveforms if gwfparams cannot be found in folder
    try
        extracted_wfs = load(fullfile(day_dir, 'extracted_wfs.mat'));
        wf = extracted_wfs.wf;
        gwfparams = extracted_wfs.gwfparams;
        clear('extracted_wfs');
        
        % Filter in good units
        if only_good
            good_clusters = gwfparams.cluster_quality.cluster_id(gwfparams.cluster_quality.group(:,1)=='g');
            [~, good_clusters_idx] = intersect(gwfparams.good_clusters, good_clusters);
            
            wf.unitIDs = wf.unitIDs(good_clusters_idx);
            wf.spikeTimeKeeps = wf.spikeTimeKeeps(good_clusters_idx, :);
            wf.waveForms = wf.waveForms(good_clusters_idx, :, :, :);
            wf.waveFormsMean = wf.waveFormsMean(good_clusters_idx, :, :);
            wf.allSpikeTimePoints = wf.allSpikeTimePoints(good_clusters_idx, :);
        end
        
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile') || strcmp(ME.identifier, 'MATLAB:nonExistentField')
            fprintf('\nextracted_wfs.mat file not found or one of the variables missing; generating a new one\n')
        else
            fprintf([ME.identifier '\n'])  % file not found has no identifier?? C'mon MatLab...
            fprintf([ME.message '\n'])
            return  % Continue here instead of break because I don't know how to catch 'file not found' exception; maybe using ME.message?
        end
    end

    if ~exist('wf', 'var') || ~exist('gwfparams', 'var')
        fprintf('Running wf extraction...\n')
        [wf, gwfparams] = get_waveforms_from_folder(day_dir, bp_filter, only_good);
        % Save gwfparams and wf for future use
        fprintf('Saving gwfparams and wf structs to mat file\n')
        save(fullfile(day_dir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
    end
end