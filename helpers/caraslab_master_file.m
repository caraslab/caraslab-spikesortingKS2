%Make sure that kilosort folder and npy-matlab folder are on MATLAB path!

% Notes: kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set to allow a fraction  of all templates 
% to span more channel groups, so that they can capture shared 
% noise across all channels. 
%
% This option is ops.criterionNoiseChannels = 0.2; 
% 
% If this number is less than 1, it will be treated as a fraction of the total number of clusters
% 
% If this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group.

%--------------------------------------------------------------------------

%Define the path to the config file. Should be in same directory as master
%file.
pathToYourConfigFile = '/Users/Melissa/Desktop/KilosortTest/'; 
run(fullfile(pathToYourConfigFile, 'caraslab_StandardConfig.m'))

%Start timer
tic; 

%Initialize GPU (will erase any existing GPU arrays)
if ops.GPU     
    gpuDevice(1); 
end

%Convert data, only for OpenEphys
if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  
end

%Preprocess data and extract spikes for initialization
[rez, DATA, uproj] = preprocessData(ops); 

%Fit templates iteratively
rez                = fitTemplates(rez, DATA, uproj);  

%Extract final spike times (overlapping extraction)
rez                = fullMPMU(rez, DATA);

% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

%Save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

%Save python results file for Phy
rezToPhy(rez, ops.root);

%Remove temporary file
delete(ops.fproc);

