% linear_regression.m

%% config
max_time_lag            = 10; %sample
window_length           = 64;
train_info_name         = 'EEG';
train_signal_name       = 'EEG_downsampled';
train_groundtruth_name  = 'EEG_ground_truth';
test_info_name          = 'test_EEG';
test_signal_name        = 'test_EEG_downsampled';
test_groundtruth_name   = 'test_EEG_ground_truth';
subject_id              = 6;

%max_time_lag = window_length;
%% main program

correct = 0;
wrong = 0;

% load training_data
eval(['train_info = ',train_info_name,';']);
% target_info = target_info{subject_id, :};
eval(['train_signal = ',train_signal_name,';']);
% target_signal = target_signal{subject_id, :};

% load ground truth
% addpath('../64_channel_Biosemi_EEG_data/Envelopes');
% for measurement_idx=1:size(target_info, 2)
%     ground_truth_channel{measurement_idx} = target_info{subject_id, measurement_idx}.trial.attended_track;
%     ground_truth{measurement_idx} = load(target_info{subject_id, measurement_idx}.trial.stimuli{ground_truth_channel{measurement_idx}},'envelope1');
% end
eval(['train_ground_truth = ',train_groundtruth_name,';']);


eval(['test_info = ',test_info_name,';']);
eval(['test_signal = ',test_signal_name,';']);
eval(['test_ground_truth = ',test_groundtruth_name,';']);


% process over the time?
%time_length = 100;%min( length(ground_truth{measurement_idx}.channel_a_downsample), ...
                      %        size(target_signal{subject_id, measurement_idx},1))-max_time_lag;
time_length = min( length(train_ground_truth{subject_id}.channel_a_downsample), ...
                   size(train_signal{subject_id, 1},1))-10*max_time_lag;

test_signal_reconstruct_a = zeros(time_length, size(test_info, 2));
test_signal_reconstruct_u = zeros(time_length, size(test_info, 2));

for subject_id = 1:12

for time_idx = 1:time_length
                      %min( length(ground_truth{measurement_idx}.channel_a_downsample), ...
                      %        size(target_signal{subject_id, measurement_idx},1))-max_time_lag;
%     decoder_a{measurement_idx} = zeros(max_time_lag*size(target_signal{subject_id, measurement_idx},2),...
%                                      max_time_lag*size(target_signal{subject_id, measurement_idx},2),...
%                                      time_length);
%     decoder_u{measurement_idx} = decoder_a{measurement_idx};
    decoder_a = zeros(max_time_lag*size(train_signal{subject_id, 1},2),...
                      1,...
                      size(train_signal, 2));
    decoder_u = decoder_a;
    for measurement_idx=1:size(train_info, 2)
        % generate lag matrix
        m = zeros(max_time_lag*size(train_signal{subject_id, measurement_idx},2), 1);
        for lag_idx = 1:max_time_lag
            for channel_idx = 1:size(train_signal{subject_id, measurement_idx},2)
                m((channel_idx-1)*max_time_lag+lag_idx) = train_signal{subject_id, measurement_idx}(time_idx+lag_idx-1, channel_idx);
            end
        end
        
        % get auto-corr and cross-corr
        m_auto_corr = m*transpose(m);
        
        % truncate result
%         audio_attend = ground_truth{subject_id, measurement_idx}.channel_a_downsample(time_idx:time_idx+max_time_lag);
%         audio_unattend = ground_truth{subject_id, measurement_idx}.channel_u_downsample(time_idx:time_idx+max_time_lag);
        
        % build except decoder result
        audio_attend = m .* train_ground_truth{subject_id, measurement_idx}.channel_a_downsample(time_idx);
        audio_unattend = m .* train_ground_truth{subject_id, measurement_idx}.channel_u_downsample(time_idx);
        
        % linear_regressive...
        decoder_a_new = lsqminnorm(audio_attend, m_auto_corr);
        decoder_u_new = lsqminnorm(audio_unattend, m_auto_corr);
        
%         decoder_a{measurement_idx}(:,:,time_idx) = decoder_a_new;
%         decoder_u{measurement_idx}(:,:,time_idx) = decoder_u_new;
        decoder_a(:,:,measurement_idx) = decoder_a_new;
        decoder_u(:,:,measurement_idx) = decoder_u_new;
    end
    
    decoder_a_avg = mean(decoder_a,3);
    decoder_u_avg = mean(decoder_u,3);
    
    for measurement_idx=1:size(test_info, 2)
        % generate lag matrix
        m = zeros(max_time_lag*size(test_signal{subject_id, measurement_idx},2), 1);
        for lag_idx = 1:max_time_lag
            for channel_idx = 1:size(test_signal{subject_id, measurement_idx},2)
                m((channel_idx-1)*max_time_lag+lag_idx) = test_signal{subject_id, measurement_idx}(time_idx+lag_idx-1, channel_idx);
            end
        end
        
        % get auto-corr and cross-corr
        m_auto_corr = m*transpose(m);

         temp = m_auto_corr*decoder_a_avg./m;
         test_signal_reconstruct_a(time_idx, measurement_idx) = temp(1);
         
         temp = m_auto_corr*decoder_u_avg./m;
         test_signal_reconstruct_u(time_idx, measurement_idx) = temp(1);
    end
end



for measurement_idx=1:size(test_info, 2)
    disp(['For data ' num2str(measurement_idx)])
    
    temp_a_1 = corrcoef(test_signal_reconstruct_a(:, measurement_idx),test_EEG_ground_truth{subject_id,measurement_idx}.channel_1_downsample(1:length(test_signal_reconstruct_a(:, measurement_idx))));
    temp_a_2 = corrcoef(test_signal_reconstruct_a(:, measurement_idx),test_EEG_ground_truth{subject_id,measurement_idx}.channel_2_downsample(1:length(test_signal_reconstruct_a(:, measurement_idx))));
    if abs(temp_a_1(1,2))>abs(temp_a_2(1,2))
        channel_a = 1;
    else
        channel_a = 2;
    end
    disp(['Test_a similar to channel ' num2str(channel_a) '(' num2str(temp_a_1(1,2)) ',' num2str(temp_a_2(1,2)) ')'])
       
    temp_u_1 = corrcoef(test_signal_reconstruct_u(:, measurement_idx),test_EEG_ground_truth{subject_id,measurement_idx}.channel_1_downsample(1:length(test_signal_reconstruct_u(:, measurement_idx))));
    temp_u_2 = corrcoef(test_signal_reconstruct_u(:, measurement_idx),test_EEG_ground_truth{subject_id,measurement_idx}.channel_2_downsample(1:length(test_signal_reconstruct_u(:, measurement_idx))));
    if abs(temp_u_1(1,2))>abs(temp_u_2(1,2))
        channel_u = 1;
    else
        channel_u = 2;
    end
    disp(['Test_u similar to channel ' num2str(channel_u) '(' num2str(temp_u_1(1,2)) ',' num2str(temp_u_2(1,2)) ')'])
    
    if channel_a ~= channel_u
        channel = channel_a;
    else
        if max([abs(temp_a_1(1,2)),abs(temp_a_2(1,2)),abs(temp_u_1(1,2)),abs(temp_u_2(1,2))])==max([abs(temp_a_1(1,2)),abs(temp_u_2(1,2))])
            channel = 1;
        else
            channel = 2;
        end
    end
    
    result{subject_id,measurement_idx}.corr = [temp_a_1(1,2),temp_a_2(1,2),temp_u_1(1,2),temp_u_2(1,2)];
    result{subject_id,measurement_idx}.channel = channel;
    
    disp(['The predict answer is attend to channel ' num2str(channel)])
    
    disp(['The correct answer is attend to channel ' num2str(test_info{subject_id, measurement_idx}.trial.attended_track)])
    
    if channel==test_info{subject_id, measurement_idx}.trial.attended_track
        correct = correct+1;
    else
        wrong = wrong+1;
    end
    
end

end

disp(['Acc:' num2str(correct/(correct+wrong)*100) '%'])
