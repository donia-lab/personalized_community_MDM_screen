function unique_element_names = find_unique_elements(samples,manifest,level)
%This function finds taxonomic elements unique to a given set of samples

sample_manifest = manifest(ismember(manifest.sample,samples),:);
sample_elements = sample_manifest{:,level};

non_sample_manifest = manifest(~ismember(manifest.sample,samples),:);
non_sample_elements = non_sample_manifest{:,level};

%Find elements present in at least one of the sample set
sample_presence = zeros(length(sample_elements{1}.Variables),1);
for i = 1:length(sample_elements)
    sample_presence = sample_presence | sample_elements{i}.Variables > 0;
end

%Find elements present in at least one of the non-sample set
non_sample_presence = zeros(length(non_sample_elements{1}.Variables),1);
for i = 1:length(non_sample_elements)
    non_sample_presence = non_sample_presence | non_sample_elements{i}.Variables > 0;
end

unique_element_index = (sample_presence) & (~non_sample_presence);
element_names = sample_elements{1}.Properties.RowNames;
unique_element_names = element_names(unique_element_index);


end

