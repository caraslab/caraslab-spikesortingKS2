function caraslab_downsample_continuous_data(input_dir,output_dir,newFs,sel)
% caraslab_downsample_continuous_data(input_dir,output_dir,newFs,sel)

if nargin < 3 || isempty(newFs), newFs = 600; end
if nargin < 4 || isempty(sel), sel = 0; end




%Check that tank directory exists and abort if it doesn't
if ~exist(input_dir,'dir')
    fprintf('\n Data directory does not exist!!\n')
    return
end



%Check if save directory exists. If it doesn't, create it now.
if ~exist(output_dir,'dir')
    [success,message,messageID] = mkdir(output_dir);
    
    %Stop if directory cannot be created, and display reason why
    if ~success
        message %#ok<*NOPRT>
        messageID
        return
    end   
end



if sel
    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(input_dir,'Select data directory');
    blocknames = {};
    for i=1:length(datafolders_names)
        [~, blocknames{end+1}, ~] = fileparts(datafolders_names{i});
    end
    
else
    %Get a list of all BLOCKS in the tank directory
    blocks = caraslab_lsdir(input_dir);
    blocknames = {blocks.name};
    
end



%Check that at least one block has been selected
if isempty(blocknames)
    fprintf('\n No BLOCKS could be found!!\n')
    return
end

cellfun(@(a) preprocesslfp(a,input_dir,output_dir,newFs),blocknames)





function preprocesslfp(block,indir,outdir,newFs)


d = dir(fullfile(indir,block,'**/settings.xml'));
pthSessionRecNode = d.folder;


sessionInfo = get_session_info(pthSessionRecNode);
 
% Find the index of the recording node (assume only one)
indNode = contains(sessionInfo.processors(:,2),'Record Node');

channelsAll = sessionInfo.processors{indNode, 3}{1};
channelsPhys = channelsAll(contains(channelsAll,'CH'));

newinfo.block = block;
newinfo.indir = indir;
newinfo.outdir = outdir;
newinfo.pthSessionRecNode = pthSessionRecNode;
newinfo.sessionInfo = sessionInfo;
newinfo.sampleRate = newFs;

data = cell(length(channelsPhys),1);
parfor i = 1:length(channelsPhys)
    ffn = fullfile(pthSessionRecNode,channelsPhys{i});
    
    fprintf('Downsampling channel %2d of %d\n',i,length(channelsPhys))
    [signal,~,info] = load_open_ephys_data_faster(ffn);
    
    origFs = info.header.sampleRate;
    
    % downsample now and filter later
    data{i} = single(resample(signal,newFs,origFs,10));
    
end

data = cell2mat(data');

outfn = fullfile(outdir,block,sprintf('%s_LFP.mat',block));

fprintf('Saving "%s" ...',outfn)

info = newinfo;
save(outfn,'data','info');

fprintf(' done\n')
























