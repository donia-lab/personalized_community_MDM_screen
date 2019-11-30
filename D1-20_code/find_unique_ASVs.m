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



%% Verify unique ASVs and see what samples they are present in

test_asv = '1c2edba6f16fd68a4972f9933d463941';
samples = filtered_manifest.sample;
test_asv_abun = zeros(length(samples),1);
for i = 1:length(samples)
    test_asv_abun(i) = filtered_manifest.rel_asv{i}{test_asv,:};
end

present_samples = samples(test_asv_abun > 0);