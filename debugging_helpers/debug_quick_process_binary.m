Tankdir = '/mnt/132bfc10-ead6-48da-986e-007a5a3d1d87/TDT tank/SUBJ-ID-150';
Savedir = '/mnt/132bfc10-ead6-48da-986e-007a5a3d1d87/TDT tank/SUBJ-ID-150';
BLOCKNAMES = {'2021-04-13_16-08-04__BasicCharacterization4'};
nChansTotal = 64;
nDataChannels = 64;
adc_channel = 0;

for i = 1:numel(BLOCKNAMES)
    
    cur_path.name = BLOCKNAMES{i};

    cur_savedir = [Savedir filesep cur_path.name];
    channels_to_extract = 1:nChansTotal;

    FULLPATH = fullfile(Tankdir,cur_path.name);


    fprintf('\n======================================================\n')
    fprintf('Processing ephys data, %s.......\n', cur_path.name)

    % Assume only one experiment folder per recording
    recording_dir = dir(fullfile(FULLPATH, '**', 'continuous.dat'));

    fidIn = fopen(fullfile(recording_dir.folder, recording_dir.name), 'r');
    fidOut_data = fopen([fullfile(recording_dir.folder, recording_dir.name) '_dataChannels.dat'], 'w');
    fidOut_adc = fopen([fullfile(recording_dir.folder, recording_dir.name) '_adcChannel.dat'], 'w');

    chunkInd = 1;
    chunkSize = 1000000;
    while 1
        cur_dat = fread(fidIn, [nChansTotal chunkSize], '*int16');
        if ~isempty(cur_dat)
            fwrite(fidOut_data, cur_dat(1:nDataChannels, :), 'int16'); % append to .dat file
            fwrite(fidOut_adc, cur_dat(adc_channel, :), 'int16'); % append to .dat file
        else
            break
        end
    end
    fclose(fidIn);

    fclose(fidOut_data);
    fclose(fidOut_adc);
end