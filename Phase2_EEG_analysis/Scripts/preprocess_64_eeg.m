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
        EEG_downsampled{subject_idx, data_idx} = transpose(resample(double(transpose(EEG{subject_idx, data_idx}.trial.RawData.EegData)), fs, EEG{subject_idx, data_idx}.trial.FileHeader.SampleRate));
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
        test_EEG_downsampled{subject_idx, data_idx} = transpose(resample(double(transpose(test_EEG{subject_idx, data_idx}.trial.RawData.EegData)), fs, test_EEG{subject_idx, data_idx}.trial.FileHeader.SampleRate));
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
