function [manifest,compound_list] = import_profinder(manifest_loc,data_loc)

%Import manifest and rename samples to match matlab convention
manifest = readtable(manifest_loc);
manifest.name = strrep(manifest.name,'-','_');

%Replace NaNs with zeros
data = readtable(data_loc);
temp_data = data.Variables;
temp_data(isnan(temp_data)) = 0;
data.Variables = temp_data;

%Generate compound list and remove compound info from data matrix
compound_list = [data.Mass,data.RT];
data = data(:,6:end);

%Get actual number of samples in the file
present_samples = intersect(manifest.name,data.Properties.VariableNames);
present_sample_num = length(present_samples); 
manifest = manifest(contains(manifest.name,present_samples),:);

%Fill manifest with compound areas
entity_num = size(compound_list,1);
entity_names = cellstr(strcat('entity_',string(1:entity_num)));
new_entities = array2table(NaN(present_sample_num,entity_num),'VariableNames',...
    entity_names);

for i = 1:present_sample_num
    sample_name = manifest.name{i}; 
    sample_index = strcmp(data.Properties.VariableNames,sample_name);
    new_entities{i,:} = transpose(data(:,sample_index).Variables);
end

manifest = [manifest, new_entities];

end

