% preprocess_24_eeg.m

%% config
eeg_file_folder = '../24_channel_Smarting_EEG_data/';

%% main process
% some static config
Number_of_subjects = 12;
Number_of_gdf_per_person = 4;


% add third_party
%addpath('../third_party/biosig4c++-2.0.2_mex');
addpath('../third_party/eeglab2019_1');

% load eeg data
for subject_idx = 1:Number_of_subjects
    current_dir = dir(strcat(eeg_file_folder,'Subject',num2str(subject_idx),'/'));
    for file_idx = 1:Number_of_gdf_per_person
        EEG{subject_idx, file_idx} = pop_biosig(strcat(eeg_file_folder,'Subject',num2str(subject_idx),'/',current_dir(file_idx+2).name));
    end
end

% eegread?????

% downsample again
% audio = resample(audio, 20, audio_fs);
% audio_fs = 20;

% BPF
% [b,a] = butter(5,[1 9]/(audio_fs/2),'bandpass');
% audio = fftfilter(b,a,audio);  %filtered signal

