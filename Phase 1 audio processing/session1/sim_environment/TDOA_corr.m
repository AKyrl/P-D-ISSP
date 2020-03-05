% TDOA_corr.m

%% config
run create_micsigs.m
run config.m
 

%% xcorr 

source{1} = original_source{1};
source{2} = original_source{2};
% source{3} = original_noise{1};
% source{4} = original_noise{2};

for j=1:size(Mic, 1)
    for i=1:length(source)
        [a, b] = xcorr(source{i}, squeeze(Mic(j, i, :)));
        TDOA{j,i} = b(find(abs(a)==max(abs(a))));
    end
end