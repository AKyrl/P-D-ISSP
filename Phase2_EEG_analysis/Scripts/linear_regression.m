% linear_regression.m

%% config
max_time_lag        = 10; %sample
window_length       = 64;
ECG_info_name       = 'EEG';
ECG_signal_name     = 'EEG_downsampled';
ECG_groundtruth_name= 'EEG_ground_truth';
subject_id          = 10;



%max_time_lag = window_length;
%% main program

% load training_data
eval(['target_info = ',ECG_info_name,';']);
% target_info = target_info{subject_id, :};
eval(['target_signal = ',ECG_signal_name,';']);
% target_signal = target_signal{subject_id, :};

% load ground truth
% addpath('../64_channel_Biosemi_EEG_data/Envelopes');
% for measurement_idx=1:size(target_info, 2)
%     ground_truth_channel{measurement_idx} = target_info{subject_id, measurement_idx}.trial.attended_track;
%     ground_truth{measurement_idx} = load(target_info{subject_id, measurement_idx}.trial.stimuli{ground_truth_channel{measurement_idx}},'envelope1');
% end
eval(['ground_truth = ',ECG_groundtruth_name,';']);


% process over the time?
for measurement_idx=1:size(target_info, 2)
    time_length = 100;%min( length(ground_truth{measurement_idx}.channel_a_downsample), ...
                      %        size(target_signal{subject_id, measurement_idx},1))-max_time_lag;
    decoder_a{measurement_idx} = zeros(max_time_lag*size(target_signal{subject_id, measurement_idx},2),...
                                     max_time_lag*size(target_signal{subject_id, measurement_idx},2),...
                                     time_length);
    decoder_u{measurement_idx} = decoder_a{measurement_idx};
    for time_idx = 1:time_length
        % generate lag matrix
        m = zeros(max_time_lag*size(target_signal{subject_id, measurement_idx},2), 1);
        for lag_idx = 1:max_time_lag
            for channel_idx = 1:size(target_signal{subject_id, measurement_idx},2)
                m((channel_idx-1)*max_time_lag+lag_idx) = target_signal{subject_id, measurement_idx}(time_idx+lag_idx-1, channel_idx);
            end
        end
        
        % get auto-corr and cross-corr
        m_auto_corr = m*transpose(m);
        
        % truncate result
%         audio_attend = ground_truth{subject_id, measurement_idx}.channel_a_downsample(time_idx:time_idx+max_time_lag);
%         audio_unattend = ground_truth{subject_id, measurement_idx}.channel_u_downsample(time_idx:time_idx+max_time_lag);
        
        % build except decoder result
        audio_attend = m_auto_corr .* ground_truth{subject_id, measurement_idx}.channel_a_downsample(time_idx);
        audio_unattend = m_auto_corr .* ground_truth{subject_id, measurement_idx}.channel_u_downsample(time_idx);
        
        % linear_regressive...
        decoder_a_new = lsqminnorm(audio_attend, m_auto_corr);
        decoder_u_new = lsqminnorm(audio_unattend, m_auto_corr);
        
        decoder_a{measurement_idx}(:,:,time_idx) = decoder_a_new;
        decoder_u{measurement_idx}(:,:,time_idx) = decoder_u_new;
    end
end




