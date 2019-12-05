%This script aggregates the results for MS1 statistical testing on all
%twenty donors
%% Aggregate and test for new metabolites

clear;clc
donor_num = 20;
data_folder = 'profinder_results/';
donor_manifest_loc = [data_folder,'20x20_donor_manifest.xlsx'];
BG_manifest_loc = [data_folder,'20x20_BG_manifest.xlsx'];
HK_manifest_loc = [data_folder,'20x20_HK_manifest.xlsx'];
data_loc_template = '_profinder.csv';

results_folder = 'untargeted_metabolomics_results/';
load_existing_analysis = true;

large_metabolite_table = [];
for i = 1:donor_num
    data_loc = [data_folder,'D',num2str(i),data_loc_template];
    putative_metabolite_table = process_untargeted_metabolomics(donor_manifest_loc,...
        BG_manifest_loc,HK_manifest_loc,data_loc,i,results_folder,load_existing_analysis);
    putative_metabolite_table.Properties.VariableNames = {['D',num2str(i)]};
    large_metabolite_table = [large_metabolite_table, putative_metabolite_table];
end


%% Remove known metabolites from HD1 screen

mass_cutoff = 0.01;
rt_cutoff = 0.2;
all_cutoff = [mass_cutoff, rt_cutoff];

novel_metabolite_table = large_metabolite_table;
novel_metabolite_table.Variables = cell(size(large_metabolite_table));
exclusion_list = readtable('targeted_compounds.csv','ReadRowNames',true);

for i = 1:size(large_metabolite_table,1)
    drug = large_metabolite_table.Properties.RowNames{i};
    drug_compound_index = find(contains(exclusion_list.Properties.RowNames,drug));
    drug_compounds = exclusion_list{drug_compound_index,{'MW','RT'}};
    
    if ~isempty(drug_compounds)
        for j = 1:size(large_metabolite_table,2)
            donor_metabolites = large_metabolite_table{drug,j}{1};
            for k = 1:size(drug_compounds,1)
                temp_metabolites = repmat(drug_compounds(k,:),size(donor_metabolites,1),1);
                full_cutoff = repmat(all_cutoff,size(donor_metabolites,1),1);
                found_metabolite_index = ...
                    abs(donor_metabolites - temp_metabolites) < full_cutoff;
                found_metabolite_index = found_metabolite_index(:,1) & found_metabolite_index(:,2);
                donor_metabolites = donor_metabolites(~found_metabolite_index,:);
                if sum(found_metabolite_index) > 0
                    exclusion_list{drug_compound_index(k),'found'} = ...
                        exclusion_list{drug_compound_index(k),'found'} + 1;
                end
            end
            novel_metabolite_table{drug,j}{1} = donor_metabolites;
        end
    end
    
end


%% Select out all metabolites and where they appear

unique_table = [];
unique_met_list = [];
var_names = {'drug','mass','RT','donors'};
base_matrix = {'string',[],[],[]};
multiple_drug_met = [];
%Loop through all metabolites in all donor/drug combinations
for i = 1:size(novel_metabolite_table,1)
    for j = 1:size(novel_metabolite_table,2)
        met_table = novel_metabolite_table{i,j}{1};
        for k = 1:size(met_table,1)
            current_met = met_table(k,:);
            full_cutoff = repmat(all_cutoff,size(unique_met_list,1),1);
            compound = novel_metabolite_table.Properties.RowNames{i};
            donor = novel_metabolite_table.Properties.VariableNames{j};
            donor = str2num(donor(2:end));
            temp_table = cell2table(base_matrix,...
                'VariableNames',var_names);
            temp_table.drug = {compound};
            temp_table.mass = current_met(1);
            temp_table.RT = current_met(2);
            temp_table.donors = {donor};
            
            %If there are currently metabolites, search against them
            if ~isempty(unique_met_list)
                matching_met = ...
                    abs(unique_met_list - repmat(current_met,size(unique_met_list,1),1)) < full_cutoff;
                matching_met = matching_met(:,1) & matching_met(:,2);
                %If the metabolites match, add its donor to the list
                if sum(matching_met) == 1
                    prior_donors = unique_table.donors(matching_met);
                    prior_donors = prior_donors{1};
                    unique_table.donors(matching_met) = {unique([prior_donors; donor])};
                    %If the metabolite's compound doesn't match the initial
                    %one, mark it for removal at the end
                    if ~strcmp(unique_table.drug{matching_met},compound)
                        multiple_drug_met = [multiple_drug_met,find(matching_met)];
                    end
                %If the metabolite doesn't match add it to the list    
                else
                    unique_met_list = [unique_met_list; current_met];
                    unique_table = [unique_table; temp_table];
                end
            %If there are no metabolites, add it to the table    
            else
                unique_met_list = [unique_met_list; current_met];
                unique_table = [unique_table; temp_table];
            end
        end
    end
end

%Remove multiply present metabolites
multiple_drug_met_index = zeros(size(unique_table,1),1);
multiple_drug_met_index(multiple_drug_met) = 1;
unique_table = unique_table(~multiple_drug_met_index,:);

%Add some other columns
unique_table.mz = unique_table.mass + 1.007;
unique_table.drug_mz = exclusion_list{unique_table.drug,'product_ion'};
unique_table.drug_RT = exclusion_list{unique_table.drug,'RT'};

%% Figure out optimal set of samples to run through molecular networking

hitting_set_table = find_optimal_sample_set(unique_table);

total_samples = cell2mat(cellfun(@(x) length(x),hitting_set_table.hitting_sets,'UniformOutput',false));

writetable(hitting_set_table,'GNPS_analysis/hitting_donor_sets.csv')

%% Add hitting sets to primary table and save

new_sets = arrayfun(@(x) num2str(x{1}),hitting_set_table.hitting_sets,'UniformOutput',false);
unique_table.donors_run = cell(size(unique_table,1),1);
for i = 1:size(unique_table,1)
    drug = unique_table.drug{i};
    current_set = new_sets{contains(hitting_set_table.drug_list,drug)};
    unique_table.donors_run{i} = current_set;
end
save('novel_untargeted_metabolomics_table.mat','unique_table');
