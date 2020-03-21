% preprocess_64_eeg.m
clear all

%% config
eeg_file_folder = '../64_channel_Biosemi_EEG_data/';

%% main process
% some static config
Number_of_train_subjects = 16;
Number_of_mat_per_train_subjects = 6;
Number_of_test_subjects = 16;
Number_of_mat_per_test_subjects = 2;

% training data

% load eeg data
for subject_idx = 1:Number_of_train_subjects
    current_dir_name = strcat(eeg_file_folder,'train/',num2str(subject_idx+1,'%02d'),'_subject/');
    current_dir = dir(current_dir_name);
    for file_idx = 1:Number_of_mat_per_train_subjects
        EEG{subject_idx, file_idx} = load(strcat(current_dir_name,current_dir(file_idx+2).name),'trial');
    end
end
% downsample
fs = 20;
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_mat_per_train_subjects
        EEG_downsampled{subject_idx, data_idx} = resample(double(EEG{subject_idx, data_idx}.trial.RawData.EegData), fs, EEG{subject_idx, data_idx}.trial.FileHeader.SampleRate);
    end
end
% BPF
[b,a] = butter(5,[1 9]/(fs/2),'bandpass');
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_mat_per_train_subjects
        for channel_idx = 1:size(EEG{subject_idx, data_idx}, 1)
            EEG_downsampled{subject_idx, data_idx}(channel_idx, :) = filter(b,a,EEG_downsampled{subject_idx, data_idx}(channel_idx, :));
        end
    end
end
% load_ground_truth
addpath('../64_channel_Biosemi_EEG_data/Envelopes');
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_mat_per_train_subjects
        EEG_ground_truth{subject_idx, data_idx}.channel_id = EEG{subject_idx, data_idx}.trial.attended_track;
        EEG_ground_truth{subject_idx, data_idx}.channel_a = load(EEG{subject_idx, data_idx}.trial.stimuli{EEG_ground_truth{subject_idx, data_idx}.channel_id},'envelope1');
        EEG_ground_truth{subject_idx, data_idx}.channel_u = load(EEG{subject_idx, data_idx}.trial.stimuli{3-EEG_ground_truth{subject_idx, data_idx}.channel_id},'envelope1');
    end
end
% downsample_ground truth
fs = 20;
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_mat_per_train_subjects
        EEG_ground_truth{subject_idx, data_idx}.channel_a_downsample = resample(double(EEG_ground_truth{subject_idx, data_idx}.channel_a.envelope1), fs, 70);
        EEG_ground_truth{subject_idx, data_idx}.channel_u_downsample = resample(double(EEG_ground_truth{subject_idx, data_idx}.channel_u.envelope1), fs, 70);
    end
end

% testing data

% load eeg data
for subject_idx = 1:Number_of_test_subjects
    current_dir_name = strcat(eeg_file_folder,'train/',num2str(subject_idx+1,'%02d'),'_subject/');
    current_dir = dir(current_dir_name);
    for file_idx = 1:Number_of_mat_per_test_subjects
        test_EEG{subject_idx, file_idx} = load(strcat(current_dir_name,current_dir(file_idx+2).name),'trial');
    end
end
% downsample
fs = 20;
for subject_idx = 1:Number_of_test_subjects
    for data_idx = 1:Number_of_mat_per_test_subjects
        test_EEG_downsampled{subject_idx, data_idx} = resample(double(test_EEG{subject_idx, data_idx}.trial.RawData.EegData), fs, test_EEG{subject_idx, data_idx}.trial.FileHeader.SampleRate);
    end
end
% BPF
[b,a] = butter(5,[1 9]/(fs/2),'bandpass');
for subject_idx = 1:Number_of_test_subjects
    for data_idx = 1:Number_of_mat_per_test_subjects
        for channel_idx = 1:size(EEG{subject_idx, data_idx}, 1)
            test_EEG_downsampled{subject_idx, data_idx}(channel_idx, :) = filter(b,a,test_EEG_downsampled{subject_idx, data_idx}(channel_idx, :));
        end
    end
end
% load ground truth
for subject_idx = 1:Number_of_test_subjects
    for data_idx = 1:Number_of_mat_per_test_subjects
        test_EEG_ground_truth{subject_idx, data_idx}.channel_id = test_EEG{subject_idx, data_idx}.trial.attended_track;
        test_EEG_ground_truth{subject_idx, data_idx}.channel_a = load(test_EEG{subject_idx, data_idx}.trial.stimuli{test_EEG_ground_truth{subject_idx, data_idx}.channel_id},'envelope1');
        test_EEG_ground_truth{subject_idx, data_idx}.channel_u = load(test_EEG{subject_idx, data_idx}.trial.stimuli{3-test_EEG_ground_truth{subject_idx, data_idx}.channel_id},'envelope1');
    end
end
% downsample_ground truth
fs = 20;
for subject_idx = 1:Number_of_test_subjects
    for data_idx = 1:Number_of_mat_per_test_subjects
        test_EEG_ground_truth{subject_idx, data_idx}.channel_a_downsample = resample(double(test_EEG_ground_truth{subject_idx, data_idx}.channel_a.envelope1), fs, 70);
        test_EEG_ground_truth{subject_idx, data_idx}.channel_u_downsample = resample(double(test_EEG_ground_truth{subject_idx, data_idx}.channel_u.envelope1), fs, 70);
    end
end