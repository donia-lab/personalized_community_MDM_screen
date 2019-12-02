% This script identifies ASVs unique to a given donor (i.e. that are unique
% to a given donor in the fecal samples and unique to that same donor in
% the cultures)

%% Import data

clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%% Go through each donor and find the unique elements

feces_manifest = filtered_manifest(strcmp(filtered_manifest.media,'feces'),:);
culture_manifest = filtered_manifest(~strcmp(filtered_manifest.media,'feces'),:);

num_unique_ASVs = zeros(max(filtered_manifest.donor),1);
unique_ASVs = cell(max(filtered_manifest.donor),1);

for i = 1:max(filtered_manifest.donor)
    donor_feces_samples = feces_manifest.sample(feces_manifest.donor == i);
    donor_culture_samples = culture_manifest.sample(culture_manifest.donor == i);
    
    unique_among_feces = find_unique_elements(donor_feces_samples,feces_manifest,'asv');
    unique_among_cultures = find_unique_elements(donor_culture_samples,culture_manifest,'asv');
    
    unique_ASVs{i} = intersect(unique_among_feces,unique_among_cultures);
    num_unique_ASVs(i) = length(unique_ASVs{i});
    
end

combined_unique_ASVs = vertcat(unique_ASVs{:});
total_num_unique_ASVs = length(combined_unique_ASVs);

mean_num_unique_ASVs = mean(num_unique_ASVs);


%% Load taxonomy and make table of unique ASVs for export 

unique_table = readtable('16S_data/taxonomy.csv');
unique_table = unique_table(contains(unique_table.FeatureID,combined_unique_ASVs),:);
unique_table.donor = nan(size(unique_table.Confidence));
unique_table.conditions = cell(size(unique_table.Confidence));

for i = 1:size(unique_table,1)
    ASV = unique_table.FeatureID{i};
    test_ASV_abun = zeros(size(filtered_manifest.donor));
    for j = 1:length(test_ASV_abun)
        test_asv_abun(j) = filtered_manifest.rel_asv{j}{ASV,:};
    end
    unique_table.conditions{i} = join(unique(filtered_manifest.media(test_asv_abun > 0)),',');
    unique_table.donor(i) = unique(filtered_manifest.donor(test_asv_abun > 0));
end

unique_table = sortrows(unique_table,'donor','ascend');
writetable(unique_table,'saved_analyses/unique_ASV_table.csv')