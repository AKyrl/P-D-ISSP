% MUSIC_wideband.m
% clear all;
% should run create_micsigs.m before run MUSIC_wideband.m

%% config
run config.m
% run create_micsigs.m
% run load_micsigs_linear_array.m
source_position_index = 1;

% normalize
% for i=1:4
%   Mic(1,i,:)  = normalize(Mic(1,i,:));
% end

%% calculate DOA of each 11pair of mics
% for source_position_index=1:6
%     disp(['source_position_index :' impulse_positions(source_position_index) ]);
% DOA_est = MUSIC_narrowband_estimation(Mic(source_position_index, [1,2,3,4],:), 5)
% DOA_est = MUSIC_wideband_estimation(Mic(source_position_index, [1,2,3,4],:), 5)
DOA_est = MUSIC_wideband_estimation(Mic(source_position_index, :,:), 5)
DOA_est2 = 90-DOA_est;

disp(['DOA_est:' num2str(DOA_est) ]);
disp(['DOA_est:' num2str(DOA_est2) ]);

% end
function [DOA_est] = MUSIC_narrowband_estimation(target_signals, mic_distance_)
    %% config
    run config.m
    
    %% load target audio source
    % eval(['target_signals = ',target_signals_name,';']);
    target_signals = reshape(target_signals, [size(target_signals,2),size(target_signals,3)]);  % 3dim to 2dim


    %% STFT
    target_signals_stft = zeros(size(target_signals,1), STFT_L, round(size(target_signals,2)/(round((100-STFT_overlap)*STFT_L/100)))-1);
    for i=1:size(target_signals,1)
        target_signals_stft(i,:,:) = stft(target_signals(i,:), ...
            'Window',hanning(STFT_L),...
            'OverlapLength',round(STFT_overlap*STFT_L/100));
    end

    %% Find max freq bin
    W_max_index =  zeros(size(target_signals_stft,1),1);
    index_list = 1:STFT_L;
    target_signals_stft_power = zeros(size(target_signals_stft,1),size(target_signals_stft,2));
    for i=1:size(target_signals_stft,1)
        for j=1:size(target_signals_stft,2)
            target_signals_stft_power(i,j) = mean(abs(target_signals_stft(i,j,:)).^2);
        end
        W_max_index(i) = index_list(find(target_signals_stft_power(i,:)==max(target_signals_stft_power(i,:)),1,'last'));
    end

    W_max_freq = (W_max_index-STFT_L/2)*(sampling_frequency/STFT_L);


    %% evalute the pseudospectrum
    angles            = [0:0.5:180];
    recorded_signal_at_w  = zeros(size(target_signals_stft,1), size(target_signals_stft,3));
    for i=1:size(target_signals_stft,1)
        recorded_signal_at_w(i,:) = target_signals_stft(i,W_max_index(i),:);
    end

    correlation_matrix = recorded_signal_at_w*recorded_signal_at_w';
    [EigenVector EigenValue ] = eig(correlation_matrix);
    EigenVector = EigenVector(:,1:size(EigenVector,2)-number_of_source_channel);

    Music_pseudospectrum = zeros(size(angles));
    Mic_position = [0:mic_distance_:mic_distance_*(size(target_signals_stft,1)-1)];

    for i=1:length(angles)
        w = W_max_freq;
        manifold_vector = exp((cosd(angles(i))*Mic_position./340./100).*w'*1i*2*pi);
        manifold_vector = reshape(manifold_vector, [size(target_signals_stft,1) 1]);
        Music_pseudospectrum(i) = 1/(manifold_vector'*EigenVector*EigenVector'*manifold_vector);
    end

    Music_pseudospectrum_sorted = findpeaks(abs(Music_pseudospectrum));
    if sum(find(Music_pseudospectrum_sorted==max(abs(Music_pseudospectrum))))==0
        Music_pseudospectrum_sorted = [max(Music_pseudospectrum) Music_pseudospectrum_sorted];
    end
    Music_pseudospectrum_sorted = sort(Music_pseudospectrum_sorted,'descend');


    for i=1:number_of_source_channel
        DOA_est(i) = angles(find(abs(Music_pseudospectrum)==Music_pseudospectrum_sorted(i),1));
    end
end

function [DOA_est] = MUSIC_wideband_estimation(target_signals, mic_distance_)
    %% config
    run config.m

    %% load target audio source
    % eval(['target_signals = ',target_signals_name,';']);
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

        Mic_position = [0:mic_distance_:mic_distance_*(size(target_signals_stft,1)-1)];

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
    Music_pseudospectrum = sum(transpose(target_signals_stft_power).*abs(Music_pseudospectrum_of_each_bins()))./sum(target_signals_stft_power);

    Music_pseudospectrum_sorted = findpeaks(Music_pseudospectrum);
    if sum(find(Music_pseudospectrum_sorted==max(Music_pseudospectrum)))==0
        Music_pseudospectrum_sorted = [max(Music_pseudospectrum) Music_pseudospectrum_sorted];
    end
    Music_pseudospectrum_sorted = sort(Music_pseudospectrum_sorted,'descend');

    for DOA_idx=1:number_of_source_channel
    %     disp(['DOA' num2str(DOA_idx) ':' num2str(angles(find(Music_pseudospectrum==Music_pseudospectrum_sorted(DOA_idx))))]);
        DOA_est(DOA_idx) = angles(find(Music_pseudospectrum==Music_pseudospectrum_sorted(DOA_idx),1));
    end
    
%     figure
%     plot(angles,Music_pseudospectrum)
end
