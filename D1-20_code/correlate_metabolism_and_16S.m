%This script computes correlations between 16S abundance and drug and
%metabolite signal

%% Load 16S data and metabolism analysis

clear;clc

manifest = wrapped_16S_import();
filtered_manifest = filter_16S_on_read_number(manifest,1e4);

metabolism_analysis = load('saved_analyses/drug_metabolism_analysis.mat');
BG_table = metabolism_analysis.BG_table;
depletion_table = metabolism_analysis.depletion_table;
mean_met_table = metabolism_analysis.mean_met_table;
sig_depletion_table = metabolism_analysis.sig_depletion_table;
sig_met_table = metabolism_analysis.sig_met_table;
met_entropy = metabolism_analysis.met_entropy;
drug_entropy = metabolism_analysis.drug_entropy;
total_entropy = [met_entropy,drug_entropy];

%% Process the drug data in a form usable for correlation analysis
num_donors = size(depletion_table,1); 
donor_order = 1:num_donors;
donor_labels = transpose(cellstr(strcat('donor_',string(donor_order))));

depletion_table = depletion_table(donor_labels,:);
mean_met_table = mean_met_table(donor_labels,:);
drug_table = [mean_met_table,depletion_table];

%Remove things with variability less than 0.5 for all drug related
%compounds
total_entropy = total_entropy(:,drug_table.Properties.VariableNames);
entropy_cutoff = 0.5;
exclusion_index = ones(size(drug_table,2),1);
for i = 1:size(drug_table,2)
    compound = strtok(drug_table.Properties.VariableNames{i},'_');
	compound_entropy = ...
        total_entropy(:,contains(total_entropy.Properties.VariableNames,compound));
    exclusion_index(i) = max(compound_entropy.Variables) > entropy_cutoff;
end

drug_table = drug_table(:,logical(exclusion_index));

%Remove remaining compounds that are always significantly depleted
drug_max_cutoff = 0.2;
drug_table = drug_table(:,transpose(max(drug_table.Variables) > drug_max_cutoff));


%% Process 16S data in a form usable for correlation analysis
BG_ind = contains(filtered_manifest.media,'BG');
BG_manifest = filtered_manifest(BG_ind,:);

%Perform cutoffs for prevalence and biomass
prevalence_cutoff = 3; 
abundance_cutoff = 0.001; %g/L
ASV_table = compute_mean_taxon_table(BG_manifest,'rel_asv',donor_order,false,prevalence_cutoff,abundance_cutoff);
genus_table = compute_mean_taxon_table(BG_manifest,'genus',donor_order,true,prevalence_cutoff,abundance_cutoff);
species_table = compute_mean_taxon_table(BG_manifest,'species',donor_order,true,prevalence_cutoff,abundance_cutoff);
family_table = compute_mean_taxon_table(BG_manifest,'family',donor_order,true,prevalence_cutoff,abundance_cutoff);


%% Compute the correlations 
corrtype = 'Spearman';

[ASV_coeff,ASV_p] = compute_table_correlations(drug_table,ASV_table,corrtype);
[species_coeff, species_p] =  compute_table_correlations(drug_table,species_table,corrtype);
[genus_coeff, genus_p] =  compute_table_correlations(drug_table,genus_table,corrtype);
[family_coeff, family_p] =  compute_table_correlations(drug_table,family_table,corrtype);

%% Correct p values and find significant correlations

significance_cutoff = 0.01;

ASV_correlation_table = find_significant_correlations(ASV_coeff,ASV_p,significance_cutoff);
species_correlation_table = find_significant_correlations(species_coeff,species_p,significance_cutoff);
genus_correlation_table = find_significant_correlations(genus_coeff,genus_p,significance_cutoff);
family_correlation_table = find_significant_correlations(family_coeff,family_p,significance_cutoff);

%% Look at correlations on the ASV level

plotting_table = ASV_correlation_table;
secondary_table = ASV_table;
for i = 1:size(plotting_table,1)
    element1 = plotting_table.item1{i};
    element2 = plotting_table.item2{i};
    plot_correlated_items(element1,element2,drug_table,secondary_table);
end

%% Look at individual correlations seen in literature

test_species = {'Bifidobacterium_s__adolescentis','Eggerthella_s__lenta'};
test_drugs = {'hydrocortisone','digoxin'};
test_mets = {'hydrocortisone_met_1','digoxin_met_1'};

nan_matrix = nan(length(test_species),6);
lit_correlation_table = array2table(nan_matrix,'RowNames',test_species,...
    'VariableNames',{'drug','met','drug_coeff','met_coeff','drug_p','met_p'});
lit_correlation_table.drug = transpose(test_drugs);
lit_correlation_table.met = transpose(test_mets);

for i = 1:length(test_species)
    species = test_species{i};
    species_index = contains(species_coeff.Properties.VariableNames,species);
    full_species = species_coeff.Properties.VariableNames{species_index};
    drug = test_drugs{i};
    met = test_mets{i};

    lit_correlation_table.drug_coeff(i) = species_coeff{drug,full_species};
    lit_correlation_table.drug_p(i) = species_p{drug,full_species};
    lit_correlation_table.met_coeff(i) = species_coeff{met,full_species};
    lit_correlation_table.met_p(i) = species_p{met,full_species};
end
