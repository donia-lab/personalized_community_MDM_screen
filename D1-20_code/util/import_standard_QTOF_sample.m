function manifest = import_standard_QTOF_sample(manifest_loc,data_loc,istd,num_met)

%This script takes a manifest and imports data for non-DMSO samples

if ~exist('num_met')
    num_met = 2;
end

%Import manifest and rename samples to match matlab convention
manifest = readtable(manifest_loc);
manifest.name = matlab.lang.makeValidName(manifest.name);

%Remove any DMSO manifest entries
manifest = manifest(~strcmp(manifest.compound,'DMSO'),:);
total_sample_num = size(manifest,1);

%Add extra columns
manifest.drug = NaN(total_sample_num,1);
manifest.istd = NaN(total_sample_num,1);
manifest.norm_drug = NaN(total_sample_num,1);

for i = 1:num_met
    met_str = num2str(i);
    manifest{:,['met_',met_str]} = NaN(total_sample_num,1);
    manifest{:,['norm_met_',met_str]} = NaN(total_sample_num,1);
end

%Import data
opts = detectImportOptions(data_loc);
opts.VariableNamesLine = 1;
opts.VariableUnitsLine = 2;
opts.DataLine = 3;
opts.RowNamesColumn = 1;
data = readtable(data_loc,opts);

%Get rid of extra row produced by matlab
data = data(:,2:end);

%Replace NaNs with zeros
temp_data = data.Variables;
temp_data(isnan(temp_data)) = 0;
data.Variables = temp_data;

%Figure out how many samples you actually have
current_sample_num = size(data,2);
max_sample = min([current_sample_num,total_sample_num]);

%Build new manifest columns from data
for i = 1:max_sample
    sample_compound = manifest.compound{i};
    sample_name = manifest.name{i};
    manifest.drug(i) = data{sample_compound,sample_name};
    manifest.istd(i) = data{istd,sample_name};
    manifest.norm_drug(i) = manifest.drug(i)/manifest.istd(i);
    
    for j = 1:num_met
        met = ['met_', num2str(j)];
        met_name = [sample_compound,'_',met];
        met_present = sum(ismember(data.Properties.RowNames,met_name));
        if met_present
            manifest{i,met} =  data{met_name,sample_name};
            manifest{i,['norm_',met]} = manifest{i,met}/manifest.istd(i);
        end
    end
end

end