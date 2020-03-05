% MUSIC_wideband.m

% should run create_micsigs.m before run MUSIC_wideband.m

%% config
% run create_micsigs.m
run config.m
% STFT_L              = 1024;
% STFT_overlap        = 50;
% target_signals_name = 'Mic(1,:,:)';
% mic_distance        = abs(m_pos(1,2)-m_pos(2,2))*100; % (cm)
% sampling_frequency  = 44100;
% number_of_source_channel = 1;
target_signals_name = 'Mic(source_position_index, [1,3],:)';
mic_distance        = 21.5; % (cm)


%% load target audio source
eval(['target_signals = ',target_signals_name,';']);
target_signals = reshape(target_signals, [size(target_signals,2),size(target_signals,3)]);  % 3dim to 2dim

%% STFT
target_signals_stft = zeros(size(target_signals,1), STFT_L, round(size(target_signals,2)/(round((100-STFT_overlap)*STFT_L/100)))-1);
for i=1:size(target_signals,1)
    target_signals_stft(i,:,:) = stft(target_signals(i,:), ...
        'Window',hanning(STFT_L),...
        'OverlapLength',round(STFT_overlap*STFT_L/100));
end


%% evalute the pseudospectrum of each frequency bins
angles            = [0:0.5:180];
Music_pseudospectrum_of_each_bins = zeros(STFT_L/2-1, length(angles));

target_signals_stft_power = zeros(size(target_signals_stft,1),STFT_L/2-1);

for frequency_bin_idx=1+STFT_L/2:STFT_L-1
    
    w = (frequency_bin_idx-STFT_L/2)*(sampling_frequency/STFT_L);
    recorded_signal_at_w  = zeros(size(target_signals_stft,1), size(target_signals_stft,3));
    for i=1:size(target_signals_stft,1)
        recorded_signal_at_w(i,:) = target_signals_stft(i,frequency_bin_idx,:);
        target_signals_stft_power(i,frequency_bin_idx-STFT_L/2) = mean(abs(recorded_signal_at_w(i,:)).^2);
    end
    
    correlation_matrix = recorded_signal_at_w*recorded_signal_at_w';
    [EigenVector EigenValue ] = eig(correlation_matrix);
    EigenVector = EigenVector(:,1:size(EigenVector,2)-number_of_source_channel);
    
    Mic_position = [0:mic_distance:mic_distance*(size(target_signals_stft,1)-1)];
    
    for angle_idx=1:length(angles)
        manifold_vector = exp((cosd(angles(angle_idx))*Mic_position./340./100).*1i.*w'*2*pi);
        manifold_vector = reshape(manifold_vector, [size(target_signals_stft,1) 1]);
        Music_pseudospectrum_of_each_bins(frequency_bin_idx-STFT_L/2, angle_idx) = 1/(manifold_vector'*EigenVector*EigenVector'*manifold_vector);
    end
    
end

% Music_pseudospectrum = geomean(abs(Music_pseudospectrum_of_each_bins));
% Music_pseudospectrum = mean(abs(Music_pseudospectrum_of_each_bins));
% weighted mean
target_signals_stft_power = sum(target_signals_stft_power);
Music_pseudospectrum = sum(transpose(target_signals_stft_power).*abs(Music_pseudospectrum_of_each_bins))./sum(target_signals_stft_power);

Music_pseudospectrum_sorted = findpeaks(Music_pseudospectrum);
if sum(find(Music_pseudospectrum_sorted==max(Music_pseudospectrum)))==0
    Music_pseudospectrum_sorted = [max(Music_pseudospectrum) Music_pseudospectrum_sorted];
end
Music_pseudospectrum_sorted = sort(Music_pseudospectrum_sorted,'descend');

for DOA_idx=1:number_of_source_channel
    disp(['DOA' num2str(DOA_idx) ':' num2str(angles(find(Music_pseudospectrum==Music_pseudospectrum_sorted(DOA_idx))))]);
    DOA_est(DOA_idx) = angles(find(Music_pseudospectrum==Music_pseudospectrum_sorted(DOA_idx),1));
end
save('DOA_est','DOA_est');
