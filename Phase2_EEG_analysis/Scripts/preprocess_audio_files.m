% preprocess_audio_files.m

%% config
audio_file_name = '../Audio_data/part1_track1_dry.wav';

%% main
% read audio file
[audio audio_fs] = audioread(audio_file_name);

% downsample
audio = resample(audio, 8000, audio_fs);
audio_fs = 8000;

% Gamma-tone filter (wait TA upload impulse response...(09-03))
%gammatoneFiltBank = gammatoneFilterBank('NumFilters', 15, 'SampleRate', audio_fs)
%audio_decomposed = gammatoneFiltBank(audio);
audio_decomposed = zeros(size(g, 1), length(audio));
for filter_idx=1:size(g, 1)
    audio_decomposed(filter_idx, :) = fftfilt(g{filter_idx}.h, audio);
end


% power-law compression
audio_decomposed = abs(audio_decomposed).^0.6;

% downsample again
audio_decomposed_fs = 20;
temp = audio_decomposed;
audio_decomposed = zeros(ceil(size(audio_decomposed).*[1,audio_decomposed_fs/audio_fs]));
for filter_idx=1:size(audio_decomposed, 1)
    audio_decomposed(filter_idx,:) = resample(temp(filter_idx,:), audio_decomposed_fs, audio_fs);
end

% linear combine?
audio_cimbined = zeros(1, size(audio_decomposed, 2));
for weight_idx=1:length(subband_weights)
    audio_cimbined = subband_weights(weight_idx).*audio_decomposed(weight_idx,:) + audio_cimbined;
end
audio_cimbined_fs = 20;

% BPF
[b,a] = butter(5,[1 9]/(audio_cimbined_fs/2),'bandpass');
audio_cimbined = filter(b,a,audio_cimbined);  %filtered signal

