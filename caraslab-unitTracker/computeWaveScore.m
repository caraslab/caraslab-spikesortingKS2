function score = computeWaveScore(wmean, channel_map, shank_channelses)
score = cell(numel(wmean)-1,1);

% Create deepcopy of wmean to store upsampled values
wmean_upsampled = deepcopy(wmean);
% Interpolate up x10
samplingRateIncrease = 10;
for iid=1:length(wmean)
    for unit1=1:length(wmean{iid})     
        % Replace values in deepcopy with upsampled holder values
        wmean_upsampled{iid}{unit1} = nan(size(wmean{iid}{unit1}) .* [samplingRateIncrease 1]);
        for chan_idx_unit1=1:size(wmean{iid}{unit1}, 2)

            % MML edit
            cur_wf_mean = wmean{iid}{unit1}(:,chan_idx_unit1);
            
            
            if all(isnan(cur_wf_mean))
                continue
            else
                newXSamplePoints = linspace(1, length(cur_wf_mean), length(cur_wf_mean) * samplingRateIncrease);
                wmean_upsampled{iid}{unit1}(:,chan_idx_unit1) = spline(1:length(cur_wf_mean), cur_wf_mean, newXSamplePoints)';
            end

            % Old code does not work
    %         wmean{iid}{unit1} = interp1(1:32,wmean{iid}{unit1}(:),1:.1:32);
        end
    end
end


% Compute score for all channels then average
% Add weights to averaging based on channel distances
for iid=1:length(wmean_upsampled)-1
    score{iid} = nan(length(wmean_upsampled{iid}),length(wmean_upsampled{iid+1}));
    for unit1=1:length(wmean_upsampled{iid})
        for unit2=1:length(wmean_upsampled{iid+1})
            weight_sum = 0;
            cor_counter = 0;
            for chan_idx_unit1=1:size(wmean_upsampled{iid}{unit1}, 2)
                w1 = wmean_upsampled{iid}{unit1}(:,chan_idx_unit1);
                w1 = w1 / norm(w1);
                for chan_idx_unit2=1:size(wmean_upsampled{iid+1}{unit2}, 2)
                    w2 = wmean_upsampled{iid+1}{unit2}(:,chan_idx_unit2);
                    w2 = w2 / norm(w2);

                    % Calculate channel distance
                    channel_unit1 = shank_channelses{iid}{unit1}(chan_idx_unit1);
                    channel_unit2 = shank_channelses{iid+1}{unit2}(chan_idx_unit2);
                    x_unit1 = channel_map.xcoords(channel_map.chanMap == channel_unit1);
                    x_unit2 = channel_map.xcoords(channel_map.chanMap == channel_unit2);
                    y_unit1 = channel_map.ycoords(channel_map.chanMap == channel_unit1);
                    y_unit2 = channel_map.ycoords(channel_map.chanMap == channel_unit2);
                    
                    % Calculate euclidean nearness (inverse-distance) and add 1 to avoid
                    % zero-division
                    % Increase this penalty factor if you're getting to many
                    % mismatching waveforms (or the inverse)
                    penalty = 0.1;
                    euc_nearness = 1/(sqrt((x_unit1 - x_unit2)^2 + (y_unit1 - y_unit2)^2)*penalty + 1);

                    weight_sum = weight_sum + euc_nearness;
                    cor_counter = cor_counter + 1;
                    cur_score = max(xcorr(w1, w2, ceil(numel(w1)/4), 'coeff')) * euc_nearness;
                    
                    score{iid}(unit1,unit2) = sum([score{iid}(unit1,unit2), cur_score], 'all','omitnan');
                end
            end
            score{iid}(unit1,unit2) = score{iid}(unit1,unit2) ./ weight_sum;  
        end
    end
    score{iid} = atanh(score{iid});
end

% Old code
% for iid=1:length(wmean_upsampled)-1
%     score{iid} = nan(length(wmean{iid}),length(wmean{iid+1}));
%     for chan_idx=1:length(wmean{iid})
%         for unit2=1:length(wmean{iid+1})
%             w1 = wmean{iid}{chan_idx};
%             w2 = wmean{iid+1}{unit2};
%             w1 = w1 / norm(w1);
%             w2 = w2 / norm(w2);
%             score{iid}(chan_idx,unit2) = atanh(max(conv(w1,w2(end:-1:1),'same')));
%             %score{iid}(unit1,unit2) = atanh(max(xcorr(wmean{iid}{unit1},wmean{iid+1}{unit2},ceil(numel(wmean{iid}{unit1})/4),'coeff')));
%         end
%     end
% end