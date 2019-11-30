function [null_mean,null_std,null_peaks] = parameterize_null_dist(manifest_loc,data_loc,istd,istd_cutoff)

%This function uses nonzero spurious peaks to compute the null distribution
%of the statistical model. This is done by searching for compound peaks
%within samples that do not contain the compound

%Import manifest and rename samples to match matlab convention
manifest = readtable(manifest_loc);
manifest.name = matlab.lang.makeValidName(manifest.name);

%Remove any DMSO manifest entries
manifest = manifest(~strcmp(manifest.compound,'DMSO'),:);
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

types = {'';'_met_1';'_met_2'};
null_peaks = []; 

for i = 1:size(manifest,1)
    sample_name = manifest.name{i};
    compound_name = manifest.compound{i};
    compound_set = [strcat(compound_name,types); istd];
    alt_name = matlab.lang.makeValidName(sample_name);
    sample_data = data(:,alt_name);
    sample_istd = sample_data(istd,:).Variables;
    off_measurements = ~contains(data.Properties.RowNames,compound_set);
    sample_data = sample_data(off_measurements,:).Variables;
    
    if sample_istd > istd_cutoff
        null_peaks = [null_peaks; sample_data./sample_istd];
    end
end

null_mean = mean(null_peaks);
null_std = std(null_peaks);

end

