% DOA_corr.m

%% config
mic_distance = 1.3; % (cm)
% mic_distance = abs(m_pos(1,2)-m_pos(2,2))*100; % (cm)
sound_speed = 340; % (m/s)
fs = 441000;

%% xcorr 

source{1} = target_signal;
% source{2} = original_source{2};

for j=1:size(Mic, 1)
    for i=1:length(source)
        [a, b] = xcorr(source{i}, Mic{j, i});
        TDOA_org{j,i} = b(find(abs(a)==max(abs(a))));
        TDOA{j,i} = b(find(abs(a)==max(abs(a))))/fs*sound_speed*100;
    end
end

%% transfer into degree

% DOA = zeros((size(TDOA,1)-1)*(size(TDOA,1)-2)/2, size(TDOA,2));
% 
% for j=1:size(TDOA, 2)
%     DOA_index = 1;
%     for i=1:size(TDOA, 1)-1
%         for k=i+1:size(TDOA, 1)
%             %         mic_distance/(TDOA{i+1,j}-TDOA{i,j})
%             DOA(DOA_index,j) = acosd((TDOA{k,j}-TDOA{i,j})/(mic_distance*(k-i)));
%             DOA_index = DOA_index+1;
%         end
%     end
% end

% DOA = zeros(size(TDOA)-[1,0]);
% 
% for j=1:size(TDOA, 2)
%     for i=1:size(TDOA, 1)-1
%         DOA(i,j) = acosd((TDOA{i+1,j}-TDOA{i,j})/mic_distance);
%     end
% end

for j=1:size(TDOA, 2)
    % DOA of two Left mic
    DOA(1,j) = acosd((TDOA{2,j}-TDOA{1,j})/1.3);
    % DOA of two right mic
    DOA(2,j) = acosd((TDOA{4,j}-TDOA{3,j})/1.3);
    % DOA of two front mic
    DOA(3,j) = acosd((TDOA{3,j}-TDOA{1,j})/21.5);
end

%% save result in file
% DOA_results = DOA(:,1);
% DOA_est(1) = mean(DOA_results(find(DOA_results==real(DOA_results))));
% DOA_results = DOA(:,2);
% DOA_est(2) = mean(DOA_results(find(DOA_results==real(DOA_results))));
% save('DOA_est','DOA_est');
