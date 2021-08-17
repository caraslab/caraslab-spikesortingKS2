tfolder = '/mnt/CL_4TB_2/temp_tanks/SUBJ-ID-150/2021-04-20_14-33-32_PassivePre/Record Node 102/';
tfile = '101_CH60.continuous';
[cur_ch_data, ~, ~] = load_open_ephys_data_faster([tfolder tfile]);


datfolder = '/mnt/CL_4TB_2/temp_tanks/SUBJ-ID-150/2021-04-13_16-08-04__BasicCharacterization4/Record Node 102/experiment1/recording1/continuous/Intan_Rec._Controller-101.0/';
datfile = 'continuous.dat';
nChansTotal = 72;
filename = [datfolder datfile];
rawsig = [];
chunkSize = 1000000;
d = dir(filename);
nSampsTotal = d.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);
% theseInds = 0;
chunkInd = 1;
fid = fopen(filename, 'r');
while 1
    dat = fread(fid, [nChansTotal chunkSize], '*int16');

    if ~isempty(dat)
        fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
        rawsig = [rawsig dat];
    else
        break
    end
    chunkInd = chunkInd+1;    
end
fclose(fid);