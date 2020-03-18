% linear_regression.m

%% config
max_time_lag        = 20; %sample
window_length       = 1024;
ECG_info_name       = 'EEG';
ECG_signal_name     = 'EEG_downsampled';
subject_id          = 10;


%% main program

% load training_data
eval(['target_info = ',ECG_info_name,';']);
% target_info = target_info{subject_id, :};
eval(['target_signal = ',ECG_signal_name,';']);
% target_signal = target_signal{subject_id, :};

% load ground truth
addpath('../64_channel_Biosemi_EEG_data/Envelopes');
for measurement_idx=1:size(target_info, 2)
    ground_truth_channel{measurement_idx} = target_info{subject_id, measurement_idx}.trial.attended_track;
    ground_truth{measurement_idx} = load(target_info{subject_id, measurement_idx}.trial.stimuli{ground_truth_channel{measurement_idx}},'envelope1');
end

% process over the time?
for measurement_idx=1:size(target_info, 2)
    for time_idx = 1:min( length(ground_truth{measurement_idx}.envelope1), ...
                          size(target_signal{subject_id, measurement_idx},1))-max_time_lag
        % generate lag matrix
        m = zeros(max_time_lag*size(target_signal{subject_id, measurement_idx},2), 1);
        for lag_idx = 1:max_time_lag
            for channel_idx = 1:size(target_signal{subject_id, measurement_idx},2)
                m((channel_idx-1)*max_time_lag+lag_idx) = target_signal{subject_id, measurement_idx}(time_idx+lag_idx-1, channel_idx);
            end
        end
        
        % get auto-corr and cross-corr
        m_auto_corr = m*transpose(m);
    end
end




