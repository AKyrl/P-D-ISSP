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


% power-law compression
audio_decomposed = abs(audio_decomposed).^0.6;

% downsample again
audio = resample(audio, 20, audio_fs);
audio_fs = 20;

% linear combine?

% BPF
[b,a] = butter(5,[1 9]/(audio_fs/2),'bandpass');
audio = fftfilter(b,a,audio);  %filtered signal

