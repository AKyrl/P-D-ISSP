% MUSIC_narrowband.m

% should run create_micsigs.m before run MUSIC_narrowband.m

%% config
STFT_L              = 1024;
STFT_overlap        = 50;
target_signals_name = 'Mic(1,:,:)+Mic(2,:,:)';
mic_distance        = abs(m_pos(1,2)-m_pos(2,2))*100; % (cm)
sampling_frequency  = 44100;
number_of_source_channel = 2;


%% load target audio source
eval(['target_signals = ',target_signals_name,';']);
target_signals = reshape(target_signals, [size(target_signals,2),size(target_signals,3)]);  % 3dim to 2dim

%% STFT
target_signals_stft = zeros(size(target_signals,1), STFT_L, round(size(target_signals,2)/(round((100-STFT_overlap)*STFT_L/100)))-1);
for i=1:size(target_signals,1)
    target_signals_stft(i,:,:) = stft(target_signals(i,:), ...
        'Window',kaiser(STFT_L,5),...
        'OverlapLength',round(STFT_overlap*STFT_L/100));
end


%% Find max freq bin
%W_max_index =  zeros(size(target_signals_stft,1),size(target_signals_stft,3));
W_max_index =  zeros(size(target_signals_stft,1),1);
index_list = 1:STFT_L;
target_signals_stft_power = zeros(size(target_signals_stft,1),size(target_signals_stft,2));
for i=1:size(target_signals_stft,1)
    %for j=1:size(target_signals_stft,3)
        %W_max_index(i,j) = index_list(find(target_signals_stft(i,:,j).^2==max(target_signals_stft(i,:,j))));
    %end
    for j=1:size(target_signals_stft,2)
        target_signals_stft_power(i,j) = mean(abs(target_signals_stft(i,j,:)).^2);
    end
    W_max_index(i) = index_list(find(target_signals_stft_power(i,:)==max(target_signals_stft_power(i,:)),1,'last'));
end

W_max_freq = (W_max_index-STFT_L/2)*(sampling_frequency/STFT_L);


%% evalute the pseudospectrum
% angles                = [-54:0.8:234];
angles_cal            = [0:0.5:180];
angles = angles_cal;
recorded_signal_at_w  = zeros(size(target_signals_stft,1), size(target_signals_stft,3));
% correlation_matrix    = zeros(size(target_signals_stft,1), size(target_signals_stft,3), size(target_signals_stft,3));
for i=1:size(target_signals_stft,1)
    recorded_signal_at_w(i,:) = target_signals_stft(i,W_max_index(i),:);
%     correlation_matrix(i,:,:) = corrcoef(recorded_signal_at_w(i,:)'*recorded_signal_at_w(i,:));
end
% 
correlation_matrix = recorded_signal_at_w*recorded_signal_at_w';
[EigenVector EigenValue ] = eig(correlation_matrix);
EigenVector = EigenVector(:,1:size(EigenVector,2)-number_of_source_channel);
% Music_pseudospectrum = zeros(size(target_signals_stft,1),length(angles));
Music_pseudospectrum = zeros(size(angles));
Mic_position = [0:mic_distance:mic_distance*(size(target_signals_stft,1)-1)];
% Mic_position = Mic_position./340./100; % cm/(m/s)/(cm/m) ==> s
for i=1:length(angles)
    w = W_max_freq;
%     -cosd(angles(i)).*1i)*([1:mic_distance:mic_distance*size(target_signals_stft,1)]
% (-cosd(angles(i)).*1i)*([0:mic_distance:mic_distance*(size(target_signals_stft,1)-1)]);
    manifold_vector = exp((cosd(angles(i))*Mic_position./340./100).*1i.*w'*2*pi);
    manifold_vector = reshape(manifold_vector, [size(target_signals_stft,1) 1]);
%     for j=1:size(target_signals_stft,1)
%     manifold_vector findpeaks= exp([0:mic_distance:mic_distance*(size(target_signals_stft,1)-1)].*(-1i).*angle(i).*w')
%         e = reshape(EigenValue(j,:),[size(EigenValue,2),1]);
    Music_pseudospectrum(i) = 1/(manifold_vector'*EigenVector*EigenVector'*manifold_vector);
%         Music_pseudospectrum(j,i) = 1/(manifold_vector'*e*e'*manifold_vector);
%     end
end

Music_pseudospectrum_sorted = findpeaks(abs(Music_pseudospectrum));

for i=1:number_of_source_channel
    disp(['DOA' num2str(i) ':' num2str(angles_cal(find(Music_pseudospectrum==Music_pseudospectrum_sorted(i))))]);
    DOA_est(i) = angles_cal(find(Music_pseudospectrum==Music_pseudospectrum_sorted(i),1));
end
save('DOA_est','DOA_est');



% correlation_matrix = zeros(size(target_signals_stft,1), size(target_signals_stft,1));
% for i=1:size(target_signals_stft,1)
%     for j=1:size(target_signals_stft,1)
% %         coef = corrcoef(recorded_signal_at_w(i,:)*recorded_signal_at_w(i,:)',recorded_signal_at_w(j,:)*recorded_signal_at_w(j,:)');
% %         correlation_matrix(i,j) = coef(1,2);
%         correlation_matrix(i,j) = recorded_signal_at_w(i,:)*recorded_signal_at_w(j,:)';
%     end
% end
% correlation_matrix =  corrcoef(recorded_signal_at_w,recorded_signal_at_w','Rows','all');
% correlation_matrix =  corrcoef(recorded_signal_at_w','Rows','all');





