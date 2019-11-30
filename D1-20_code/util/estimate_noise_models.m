function noise_models = estimate_noise_models(donor_manifest_loc,donor_loc,istd,plt)

%This script uses data from all donors to parameterize a distribution of
%background noise and a model of measurement noise

%% Load data

donor_manifest = import_standard_QTOF_sample(donor_manifest_loc,donor_loc,istd);

donor_list = unique(donor_manifest.donor);
compound_list = unique(donor_manifest.compound);
compound_list = compound_list(~cellfun('isempty',compound_list));

%% Remove outliers based on istd amount

filtered_donor = filter_outliers_from_istd(donor_manifest);

%% Estimate parameters of null distribution

istd_cutoff = 1e6;
[null_mean,null_std] = ...
    parameterize_null_dist(donor_manifest_loc,donor_loc,istd,istd_cutoff);

%% Loop through samples and get mean-stdev pairs

mean_list = [];
stdev_list = [];

mol_types = {'norm_drug','norm_met_1','norm_met_2'};
for i = 1:length(donor_list)
    donor = donor_list(i);
    for j = 1:length(compound_list)
        compound = compound_list{j};
        samples_index = (filtered_donor.donor == donor) & ...
            strcmp(filtered_donor.compound,compound);
        samples = filtered_donor(samples_index,:);
        for k = 1:length(mol_types)
            mean_list = [mean_list; mean(samples(:,mol_types{k}).Variables)];
            stdev_list = [stdev_list; std(samples(:,mol_types{k}).Variables)];
        end
    end
end

%% Use null distribution to select likely true measurements and then model
%the relationship between mean and variance

mean_cutoff = null_mean + 5*null_std;
likely_true_measurement = ~isnan(stdev_list) & (mean_list > mean_cutoff);
filt_mean_list = mean_list(likely_true_measurement);
filt_std_list = stdev_list(likely_true_measurement);

x = [zeros(size(filt_mean_list))+1,log10(filt_mean_list)];
y = log10(filt_std_list);
coefficient = x\y;

if plt
    figure
    loglog(filt_mean_list,filt_std_list,'ko')
    hold on
    loglog(filt_mean_list,10^coefficient(1).*filt_mean_list.^(coefficient(2)),'k-')
    xlabel('Mean signal')
    ylabel('Signal standard deviation')
end

noise_models.null_mean = null_mean;
noise_models.null_std = null_std;
noise_models.a = 10^coefficient(1);
noise_models.b = coefficient(2);

end