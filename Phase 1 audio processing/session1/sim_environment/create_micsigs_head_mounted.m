% create_micsigs_head_mounted.m

% Mic{1,i} ==> L1
% Mic{2,i} ==> L2
% Mic{3,i} ==> R1
% Mic{4,i} ==> R2

% clc; clear all; close all
%% config
run config.m
% length_limit = 10; % (s)
% upsampling_FS = 441000;

%% config of impulse response
impulse_positions   = ["s30";"s60";"s90";"s-30";"s-60";"s-90"];
impulse_path        = '../Head_mounted_IRs/Head_mounted_IRs/';
L1_impulse_filename = '/HMIR_L1.wav';
L2_impulse_filename = '/HMIR_L2.wav';
R1_impulse_filename = '/HMIR_R1.wav';
R2_impulse_filename = '/HMIR_R2.wav';

%% define filenames of waveforms and load waveforms

% speech_filename{1}      = 'audio_files/speech1.wav';
% speech_filename{2}      = 'audio_files/speech2.wav';
% 
% whitenoise_filename{1}  = 'audio_files/whitenoise_signal_1.wav';
% whitenoise_filename{2}  = 'audio_files/whitenoise_signal_2.wav';

dry_filename{1}         = 'audio_files/part1_track1_dry.wav';
% dry_filename{2}         = 'audio_files/part1_track2_dry.wav';
% 
% bubble_filename         = 'audio_files/Babble_noise1.wav';

% noise_filename = dry_filename;
% noise_filename = whitenoise_filename;
% 
% [source{1}, source_Fs{1}]   = audioread(speech_filename{1});
% [source{2}, source_Fs{2}]   = audioread(speech_filename{2});
% 
% 
% [noise{1}, noise_Fs{1}]   = audioread(noise_filename{1});
% [noise{2}, noise_Fs{2}]   = audioread(noise_filename{2});

[target_signal, target_Fs]   = audioread(dry_filename{1});

%% load impulse response
for i=1:length(impulse_positions)
    [impulse_L1{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), L1_impulse_filename));
%     impulse_L1{i}  = resample(impulse_L1{i}, upsampling_FS, Fs);
    [impulse_L2{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), L2_impulse_filename));
%     impulse_L2{i}  = resample(impulse_L2{i}, upsampling_FS, Fs);
    [impulse_R1{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), R1_impulse_filename));
%     impulse_R1{i}  = resample(impulse_R1{i}, upsampling_FS, Fs);
    [impulse_R2{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), R2_impulse_filename));
%     impulse_R2{i}  = resample(impulse_R2{i}, upsampling_FS, Fs);
end


%% preprocess for waveform
% target_signal  = resample(target_signal, upsampling_FS, target_Fs);
% target_signal  = target_signal(1:length_limit* upsampling_FS);
target_signal  = target_signal(1:length_limit* Fs);

%% filter data
Mic = zeros(    length(impulse_positions), ...
                4, ...
                length_limit*fs_RIR);

for i=1:length(impulse_positions)
    Mic(i,1,:) = fftfilt(impulse_L1{i}, target_signal);
    Mic(i,2,:) = fftfilt(impulse_L2{i}, target_signal);
    Mic(i,3,:) = fftfilt(impulse_R1{i}, target_signal);
    Mic(i,4,:) = fftfilt(impulse_R2{i}, target_signal);
end

%% combine signal from each speaker
sound = zeros(  4, length_limit*fs_RIR);

for i=1:4
    sound(i,:) = squeeze(Mic(1,i,:));
end

for i=2:length(impulse_positions)
    for j=1:4
        sound(j,:) = sound(j,:) + transpose(squeeze(Mic(i,j,:)));
    end
end


