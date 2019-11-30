function [optimal_media, optimal_value] = find_highest_medium(manifest,metric)

%This script finds the average optimal media relative to a specified metric

media_list = unique(manifest.media);
media_average = zeros(size(media_list));
for i = 1:length(media_list)
    
    media = media_list{i};
    partial_manifest = manifest(contains(manifest.media,media),:);
    media_average(i) = mean(partial_manifest{:,metric});
    
end

if strcmp(metric,'JS')
    [optimal_value,optimal_index] = min(media_average);
else
    [optimal_value,optimal_index] = max(media_average);
end

optimal_media = media_list{optimal_index};

end