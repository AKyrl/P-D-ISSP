%% config.m
% Essential configuration
load Computed_RIRs.mat
length_limit                = 10; % (s)
STFT_L                      = 1024;
STFT_overlap                = 50;
target_signals_name         = 'Mic(1,:,:)';
mic_distance                = abs(m_pos(1,2)-m_pos(2,2))*100; % (cm)
sampling_frequency          = 44100;
number_of_source_channel    = 1;
LMS_step_size               = 0.1;
LMS_precision_error         = 1e-10;
fs                          = sampling_frequency;