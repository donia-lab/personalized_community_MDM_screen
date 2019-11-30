function [filtered_manifest] = filter_outliers_from_istd(manifest)
%This function takes in a manifest and filters by plate for outliers based
%on internal standard and also removes nans

plate_list = unique(manifest.plate);
filtered_manifest = [];

for i = 1:length(plate_list)
    
    %Get data for each donor
    plate_data = manifest(manifest.plate == plate_list(i),:);
    plate_istd = plate_data.istd;
    
    %Compute non-outliers based on quantiles
    plate_iqr = iqr(plate_istd);
    high_lim = quantile(plate_istd,0.75) + 3*plate_iqr;
    low_lim = quantile(plate_istd,0.25) - 3*plate_iqr;
    non_outliers = (plate_istd > low_lim) & (plate_istd < high_lim);
    %non_outliers = (plate_istd > 1e6);
    filtered_manifest = [filtered_manifest; plate_data(non_outliers,:)];
    
end

%Find rows where all measurements are nan
measurement_col = contains(filtered_manifest.Properties.VariableNames,{'drug','istd','met'});
non_nan_rows = sum(isnan(filtered_manifest(:,measurement_col).Variables),2) < sum(measurement_col);
filtered_manifest = filtered_manifest(non_nan_rows,:);

end

