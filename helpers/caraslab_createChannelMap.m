function caraslab_createChannelMap(savedir,probetype, recording_format)
%caraslab_createChannelMap(savedir,probetype)
%
%Create a channel map file for a specific probe type
%
% Input variables:
%   savedir:    directory to save channelmap file (-mat)
%  
%   Probetype:  specifies the probe style used for recordings
%   Only the following have been set up:
%               'NNBuz5x1264':      Neuronexus Buzsaki 5x12 H64LP
%               'NN4x16Poly64':     Neuronexus A4x16 Poly2 H64LP (only for Synapse)
%               'NNA4x16Lin64':     Neuronexus A4x16 Poly2 Linear H64LP
%               'NNA2x32':          Neuronexus A2x32-5mm-25-200-177 (only for Synapse)
%               'NNoptrodeLin4':    Neuronexus Qtrode-Linear

%
% Kilosort2 note: kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is ops.criterionNoiseChannels = 0.2; 
% If this number is less than 1, it will be treated as a fraction of the total number of clusters
% If this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group.



%How many channels are on your probe?
switch probetype
    case {'NNBuz5x1264','NN4x16Poly64','NNA4x16Lin64', 'NNA2x32'}
        Nchannels = 64;

    case 'NNoptrodeLin4'
        Nchannels = 4;
end


%Create the channel map (and a version that's indexed starting at zero)
chanMap = 1:Nchannels;
chanMap0ind = chanMap - 1;


switch probetype
    case 'NNBuz5x1264'
        switch recording_format
            % Synapse channel mapping
            case 'synapse'
                %Define the channel groups:
                shank1L = [61,60,63,62,58,59]; %left side of shank 1
                shank1R = [54,55,52,53,57,56]; %right side of shank 1

                shank2L = [49,48,51,50,18,17]; %left side of shank 2
                shank2R = [21,24,26,25,22,20]; %right side of shank 2

                shank3L = [32,29,30,28,31,27]; %left side of shank 3
                shank3R = [2,3,4,6,1,5];       %right side of shank 3

                shank4L = [11,10,8,7,12,14];   %left side of shank 4
                shank4R = [47,46,45,44,16,15]; %right side of shank 4

                shank5L = [40,41,42,43,39,38]; %left side of shank 5
                shank5R = [35,34,64,33,36,37]; %right side of shank 5

                extrasite1 = 23;     %extra sites on shank 3
                extrasite2 = 9;     %extra sites on shank 3
                extrasite3 = 19;     %extra sites on shank 3
                extrasite4 = 13;     %extra sites on shank 3
                
            % Intan channel mapping
            case  'intan'
                %Define the channel groups:
                shank1L = [19, 20, 17, 18, 22, 21]; %left side of shank 1
                shank1R = [26, 25, 28, 27, 23, 24]; %right side of shank 1

                shank2L = [31, 32, 29, 30, 1, 2]; %left side of shank 2
                shank2R = [6, 7, 9, 10, 5, 3]; %right side of shank 2

                shank3L = [15, 14, 13, 11, 16, 12]; %left side of shank 3
                shank3R = [49, 52, 51, 53, 50, 54];       %right side of shank 3

                shank4L = [60, 57, 55, 56, 59, 61];   %left side of shank 4
                shank4R = [33, 34, 35, 36, 63, 64]; %right side of shank 4

                shank5L = [40, 39, 38, 37, 41, 42]; %left side of shank 5
                shank5R = [45, 46, 47, 48, 44, 43]; %right side of shank 5

                extrasite1 = 8;     %extra sites on shank 3
                extrasite2 = 58;     %extra sites on shank 3
                extrasite3 = 4;     %extra sites on shank 3
                extrasite4 = 62;     %extra sites on shank 3
        end
                % This will be used below to set kcoords appropriately
                kcoords_map = {[shank1L,shank1R]; [shank2L,shank2R]; ...
                    [shank3L,shank3R]; [shank4L,shank4R]; [shank5L,shank5R]; ...
                    extrasite1; extrasite2; extrasite3; extrasite4};


                %-----------------------------------------------------------------------
                %Define the x coordinates for each channel group (in relative microns)
                %-----------------------------------------------------------------------
                %On each shank, sites are spaced in two columns, set 20 um apart
                xL = zeros(1,6);
                xR = zeros(1,6)+20;

                %Shank 1 will be defined as starting at position 0
                Xshank1L = xL;
                Xshank1R = xR;

                %Shank 2 is 200 um away from shank 1
                Xshank2L = 200+xL;
                Xshank2R = 200+xR;

                %Shank 3 is 400 um away from shank 1
                Xshank3L = 400+xL;
                Xshank3R = 400+xR;

                %Shank 4 is 600 um away from shank 1
                Xshank4L = 600+xL;
                Xshank4R = 600+xR;

                %Shank 5 is 800 um away from shank 1
                Xshank5L = 800+xL;
                Xshank5R = 800+xR;

                %The extra sites are centered on shank 3 (i.e. 10 um offset from the columns)
                Xextra1 = 10 + Xshank3L(1);
                Xextra2 = 10 + Xshank3L(2);
                Xextra3 = 10 + Xshank3L(3);
                Xextra4 = 10 + Xshank3L(4);

                % This will be used below to set xcoords appropriately
                xcoords_map = {[Xshank1L,Xshank1R]; [Xshank2L,Xshank2R]; ...
                    [Xshank3L,Xshank3R]; [Xshank4L,Xshank4R]; [Xshank5L,Xshank5R]; ...
                    Xextra1; Xextra2; Xextra3; Xextra4};


                %-----------------------------------------------------------------------
                %Define the y coordinates for each channel group (in relative microns)
                %-----------------------------------------------------------------------
                %Left column starts 55 um above the tip,
                %and extends upwards in 20 um spacing
                yL = fliplr(55:20:(55+20*5));

                %Right column starts 35 um above the tip,
                %and extends upward in 20 um spacing
                yR = fliplr(35:20:(35+20*5));


                %Shank 1
                Yshank1L = yL;
                Yshank1R = yR;

                %Shank 2
                Yshank2L = yL;
                Yshank2R = yR;

                %Shank 3
                Yshank3L = yL;
                Yshank3R = yR;

                %Shank 4
                Yshank4L = yL;
                Yshank4R = yR;

                %Shank 5
                Yshank5L = yL;
                Yshank5R = yR;

                %Extra sites start 200 um above the top most site on the left,
                %and extend upwards in 200 um spacing
                Yextra1 = yL(end)+200;
                Yextra2 = yL(end)+400;
                Yextra3 = yL(end)+600;
                Yextra4 = yL(end)+800;


                % This will be used below to set ycoords appropriately
                ycoords_map = {[Yshank1L,Yshank1R]; [Yshank2L,Yshank2R]; ...
                    [Yshank3L,Yshank3R]; [Yshank4L,Yshank4R]; [Yshank5L,Yshank5R]; ...
                    Yextra1; Yextra2; Yextra3; Yextra4};

    case 'NN4x16Poly64'
        switch recording_format
            case 'synapse'
                %Define the channel groups:
                shank1L = [59,58,61,60,63,62,56,57]; %left side of shank 1
                shank1R = [52,53,50,51,48,49,55,54]; %right side of shank 1

                shank2L = [24,21,22,20,17,18,26,25]; %left side of shank 2
                shank2R = [29,32,31,27,23,19,30,28]; %right side of shank 2

                shank3L = [3,2,1,5,9,13,4,6]; %left side of shank 3
                shank3R = [10,11,12,14,15,16,8,7];       %right side of shank 3

                shank4L = [42,43,44,45,46,47,41,40];   %left side of shank 4
                shank4R = [37,36,35,34,64,33,38,39]; %right side of shank 4
            case  'intan'
                % TODO
                fprintf('This probe type has not been set up yet for the Intan system')
        end

                % This will be used below to set kcoords appropriately
                kcoords_map = {[shank1L,shank1R]; [shank2L,shank2R]; ...
                    [shank3L,shank3R]; [shank4L,shank4R]};


                %-----------------------------------------------------------------------
                %Define the x coordinates for each channel group (in relative microns)
                %-----------------------------------------------------------------------
                %On each shank, sites are spaced in two columns, set 30 um apart
                xL = zeros(1,8);
                xR = zeros(1,8)+30;

                %Shank 1 will be defined as starting at position 0
                Xshank1L = xL;
                Xshank1R = xR;

                %Shank 2 is 200 um away from shank 1
                Xshank2L = 200+xL;
                Xshank2R = 200+xR;

                %Shank 3 is 400 um away from shank 1
                Xshank3L = 400+xL;
                Xshank3R = 400+xR;

                %Shank 4 is 600 um away from shank 1
                Xshank4L = 600+xL;
                Xshank4R = 600+xR;


                % This will be used below to set xcoords appropriately
                xcoords_map = {[Xshank1L,Xshank1R]; [Xshank2L,Xshank2R]; ...
                    [Xshank3L,Xshank3R]; [Xshank4L,Xshank4R]};

                %-----------------------------------------------------------------------
                %Define the y coordinates for each channel group (in relative microns)
                %-----------------------------------------------------------------------
                %Left column starts 73 um above the tip,
                %and extends upwards in 23 um spacing; 7 spaces between channels
                yL = fliplr(73:23:(73+23*7));

                %Right column starts 50 um above the tip,
                %and extends upward in 23 um spacing; 7 spaces between channels
                yR = fliplr(50:23:(50+23*7));

                %Shank 1
                Yshank1L = yL;
                Yshank1R = yR;

                %Shank 2
                Yshank2L = yL;
                Yshank2R = yR;

                %Shank 3
                Yshank3L = yL;
                Yshank3R = yR;

                %Shank 4
                Yshank4L = yL;
                Yshank4R = yR;


                % This will be used below to set ycoords appropriately
                ycoords_map = {[Yshank1L,Yshank1R]; [Yshank2L,Yshank2R]; ...
                    [Yshank3L,Yshank3R]; [Yshank4L,Yshank4R]};
                
    case 'NNA2x32'
        switch recording_format
            case 'synapse'
                %Define the channel groups:
                shank1 = [18, 49, 17, 48, 20, 51, 22, 50, 21, 53, 24, 52, 26, 55, ...
                    25, 54, 28, 57, 30, 56, 29, 59, 32, 58, 31, 61, 27, 60, 23, 63, 19, 62];
                shank2 = setdiff(1:64, shank1);  % fill with whatever is left out of 64 channels
            case 'intan'
                % TODO:
                fprintf('This probe type has not been set up yet for the Intan system')
        end
                % This will be used below to set kcoords appropriately
                kcoords_map = {[shank1]; [shank2]};


                %-----------------------------------------------------------------------
                %Define the x coordinates for each channel group (in relative microns)
                %-----------------------------------------------------------------------
                % Shanks are 200 um apart
                Xshank1 = zeros(1,32);
                Xshank2 = zeros(1,32)+200;

                % This will be used below to set xcoords appropriately
                xcoords_map = {[Xshank1], [Xshank2]};

                %-----------------------------------------------------------------------
                %Define the y coordinates for each channel group (in relative microns)
                %-----------------------------------------------------------------------
                % Electrodes start 50 um above the tip,
                %and extend upwards in 25 um spacing; 31 spaces between channels
                Yshank1 = fliplr(50:25:(50+25*31));
                Yshank2 = Yshank1;
                % This will be used below to set ycoords appropriately
                ycoords_map = {[Yshank1], [Yshank2]};
        
        
    case 'NNA4x16Lin64'
        switch recording_format
            case 'synapse'
                %Define the channel groups:
                shank1L = [59, 58, 61, 60, 63, 62, 56]; %left side of shank 1
                shank1R = [52, 53, 50, 51, 48, 49, 55]; %right side of shank 1

                shank2L = [24, 21, 22, 20, 17, 18, 26]; %left side of shank 2
                shank2R = [29, 32, 31, 27, 23, 19, 30]; %right side of shank 2

                shank3L = [3, 2, 1, 5, 9, 13, 4]; %left side of shank 3
                shank3R = [10, 11, 12, 14, 15, 16, 8];       %right side of shank 3

                shank4L = [42, 43, 44, 45, 46, 47, 41];   %left side of shank 4
                shank4R = [37, 36, 35, 34, 64, 33, 38]; %right side of shank 4
                extrasite1shank1 = 54;     %extra sites on shank 1
                extrasite2shank1 = 57;     %extra sites on shank 1
                extrasite3shank2 = 28;     %extra sites on shank 2
                extrasite4shank2 = 25;     %extra sites on shank 2
                extrasite5shank3 = 7;     %extra sites on shank 3
                extrasite6shank3 = 6;     %extra sites on shank 3
                extrasite7shank4 = 39;     %extra sites on shank 4
                extrasite8shank4 = 40;     %extra sites on shank 4
            case 'intan'
                %Define the channel groups:
                shank1L = [21, 22, 19, 20, 17, 18, 24]; %left side of shank 1
                shank1R = [28, 27, 30, 29, 32, 31, 25]; %right side of shank 1

                shank2L = [7, 6, 5, 3, 2, 1, 9]; %left side of shank 2
                shank2R = [14, 15, 16, 12, 8, 4, 13]; %right side of shank 2

                shank3L = [52, 49, 50, 54, 58, 62, 51]; %left side of shank 3
                shank3R = [57, 60, 59, 61, 64, 63, 55];       %right side of shank 3

                shank4L = [38, 37, 36, 35, 34, 33, 39];   %left side of shank 4
                shank4R = [43, 44, 45, 46, 47, 48, 42]; %right side of shank 4
                extrasite1shank1 = 23;     %extra sites on shank 1
                extrasite2shank1 = 26;     %extra sites on shank 1
                extrasite3shank2 = 10;     %extra sites on shank 2
                extrasite4shank2 = 11;     %extra sites on shank 2
                extrasite5shank3 = 53;     %extra sites on shank 3
                extrasite6shank3 = 56;     %extra sites on shank 3
                extrasite7shank4 = 40;     %extra sites on shank 4
                extrasite8shank4 = 41;     %extra sites on shank 4
        end
            % This will be used below to set kcoords appropriately
            kcoords_map = {[shank1L,shank1R]; [shank2L,shank2R]; ...
                [shank3L,shank3R]; [shank4L,shank4R]; ...
                extrasite1shank1; extrasite2shank1; extrasite3shank2; extrasite4shank2;
                extrasite5shank3; extrasite6shank3; extrasite7shank4; extrasite8shank4};


            %-----------------------------------------------------------------------
            %Define the x coordinates for each channel group (in relative microns)
            %-----------------------------------------------------------------------
            %On each shank, sites are spaced in two columns, set 17.32 um apart
            xL = zeros(1,7);
            xR = zeros(1,7)+17.32;

            %Shank 1 will be defined as starting at position 0
            Xshank1L = xL;
            Xshank1R = xR;

            %Shank 2 is 150 um away from shank 1
            Xshank2L = 150+xL;
            Xshank2R = 150+xR;

            %Shank 3 is 300 um away from shank 1
            Xshank3L = 300+xL;
            Xshank3R = 300+xR;

            %Shank 4 is 450 um away from shank 1
            Xshank4L = 450+xL;
            Xshank4R = 450+xR;


            % Extra sites are centered on each shank, 17.32/2 from L electrodes
            extrasite1shank1 = Xshank1L(1)+17.32/2;     %extra sites on shank 1
            extrasite2shank1 = Xshank1L(1)+17.32/2;     %extra sites on shank 1
            extrasite3shank2 = Xshank2L(1)+17.32/2;     %extra sites on shank 2
            extrasite4shank2 = Xshank2L(1)+17.32/2;     %extra sites on shank 2
            extrasite5shank3 = Xshank3L(1)+17.32/2;     %extra sites on shank 3
            extrasite6shank3 = Xshank3L(1)+17.32/2;     %extra sites on shank 3
            extrasite7shank4 = Xshank4L(1)+17.32/2;     %extra sites on shank 4
            extrasite8shank4 = Xshank4L(1)+17.32/2;     %extra sites on shank 4

            % This will be used below to set xcoords appropriately
            xcoords_map = {[Xshank1L,Xshank1R]; [Xshank2L,Xshank2R]; ...
                [Xshank3L,Xshank3R]; [Xshank4L,Xshank4R]; ...
                extrasite1shank1; extrasite2shank1; extrasite3shank2; extrasite4shank2;
                extrasite5shank3; extrasite6shank3; extrasite7shank4; extrasite8shank4};


            %-----------------------------------------------------------------------
            %Define the y coordinates for each channel group (in relative microns)
            %-----------------------------------------------------------------------
            %Left column bottom channel is 50 um above the tip,
            %and extends upward in 20 um spacing
            yL = fliplr(50:20:(50+20*6));

            %Right column starts 40 um above the tip,
            %and extends upward in 20 um spacing
            yR = fliplr(40:20:(40+20*6));

            %Shank 1
            Yshank1L = yL;
            Yshank1R = yR;

            %Shank 2
            Yshank2L = yL;
            Yshank2R = yR;

            %Shank 3
            Yshank3L = yL;
            Yshank3R = yR;

            %Shank 4
            Yshank4L = yL;
            Yshank4R = yR;

            % Extra sites on shank 1 are centered on each shank, 100 and 200 away from top
            % L site; shank 2 starts at same Y as shank 1 and the next moves up
            % 100 um; same logic for 3 and 4
            extrasite1shank1 = yL(1)+100;     %extra sites on shank 1
            extrasite2shank1 = yL(1)+200;     %extra sites on shank 1
            extrasite3shank2 = yL(1)+200;     %extra sites on shank 2
            extrasite4shank2 = yL(1)+300;     %extra sites on shank 2
            extrasite5shank3 = yL(1)+300;     %extra sites on shank 3
            extrasite6shank3 = yL(1)+400;     %extra sites on shank 3
            extrasite7shank4 = yL(1)+400;     %extra sites on shank 4
            extrasite8shank4 = yL(1)+500;     %extra sites on shank 4



            % This will be used below to set ycoords appropriately
            ycoords_map = {[Yshank1L,Yshank1R]; [Yshank2L,Yshank2R]; ...
                [Yshank3L,Yshank3R]; [Yshank4L,Yshank4R]; ...
                extrasite1shank1; extrasite2shank1; extrasite3shank2; extrasite4shank2;
                extrasite5shank3; extrasite6shank3; extrasite7shank4; extrasite8shank4};
    case 'NNoptrodeLin4'
        switch recording_format
            case 'synapse'
                % From top (closest to optical fiber) to bottom
                fprintf('This probe type has not been set up yet for the TDT system')
            case 'intan'
                firstch = 2;
                secondch = 1;
                thirdch = 3;
                fourthch = 4; 
        end

            % This will be used below to set kcoords appropriately
            kcoords_map = {firstch, secondch, thirdch, fourthch};


            %-----------------------------------------------------------------------
            %Define the x coordinates for each channel group (in relative microns)
            %-----------------------------------------------------------------------

            % This will be used below to set xcoords appropriately
            xcoords_map = {0; 0; 0; 0};


            %-----------------------------------------------------------------------
            %Define the y coordinates for each channel group (in relative microns)
            %-----------------------------------------------------------------------
            % Bottom channel is 100 um above the tip and other channels
            % exten upwards in 50 um spaces (4 ch total)
            % This will be used below to set ycoords appropriately
            ycoords_map = {250; 200; 150; 100};
    otherwise
        fprintf('\nProbe dimensions not specified!\nEdit caraslab_CreateChannelMapFile.m to add dimensions before map can be generated.\n')
        return
end


%-----------------------------------------------------------------------
%Define the k coordinates for each channel group
%-----------------------------------------------------------------------
%The k coordinates indicate the group that each channel belongs to. Nearby
%sites on a single shank might pickup activity from the same neuron, for
%instance, and thus belong to the same group, but sites that are spaced far
%apart, or on different shanks, could not possibly pick up the same unit,
%and thus should be identified as being members of different groups.
%Specifying the groups will help Kilosort's algorithm discard noisy
%templates that are shared across groups.
%
% Channels need to be mapped according to the .dat channel stream which is 1:Nchan
kcoords =  ones(Nchannels, 1);  %  Placeholders
xcoords =  ones(Nchannels, 1);
ycoords =  ones(Nchannels, 1);
for x=1:length(kcoords_map)  % loop through shanks
    kcoords(kcoords_map{x}) = x;  % Change kcoords to shank index
    xcoords(kcoords_map{x}) = xcoords_map{x};  % Map xcoords
    ycoords(kcoords_map{x}) = ycoords_map{x};  % Map ycoords
end

% Identify dead (or disconnected) channels
connected = true(Nchannels, 1); % a 'connected' channel is one that is active (not dead)

%% Save
filename = fullfile(savedir,[probetype, '_', recording_format, '.mat']);
save(filename,'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind')
fprintf('Saved channel map file: %s \n',filename);

 