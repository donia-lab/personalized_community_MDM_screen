function manifest = import_DMSO_QTOF_sample(manifest_loc,data_loc,istd)

%This script takes a manifest and imports data for DMSO samples

%Import manifest and rename samples to match matlab convention
manifest = readtable(manifest_loc);
manifest.name = matlab.lang.makeValidName(manifest.name);

%Exclude all but DMSO manifest entries
manifest = manifest(strcmp(manifest.compound,'DMSO'),:);
total_sample_num = size(manifest,1);


%Import data
opts = detectImportOptions(data_loc);
opts.VariableNamesLine = 1;
opts.VariableUnitsLine = 2;
opts.RowNamesColumn = 1;
data = readtable(data_loc,opts);

%Get rid of extra row produced by matlab
data = data(:,2:end);

%Replace NaNs with zeros
temp_data = data.Variables;
temp_data(isnan(temp_data)) = 0;
data.Variables = temp_data;

%Add extra columns
manifest.istd = NaN(total_sample_num,1);
metabolite_list = ...
    data.Properties.RowNames(contains(data.Properties.RowNames,'_met_'));

manifest{:,metabolite_list} = NaN(total_sample_num,length(metabolite_list));
manifest{:,strcat('norm_',metabolite_list)} = NaN(total_sample_num,length(metabolite_list));

%Build new manifest columns from data
for i = 1:total_sample_num
    sample_name = manifest.name{i};
    manifest.istd(i) = data{istd,sample_name};
    for j = 1:length(metabolite_list)
        met_name = metabolite_list{j};
        manifest{i,met_name} =  data{met_name,sample_name};
        manifest{i,['norm_',met_name]} = manifest{i,met_name}/manifest.istd(i);
    end
end

end