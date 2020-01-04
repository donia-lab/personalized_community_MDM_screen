% This script takes in results from GNPS and check whether the parent
% molecule and the putative metabolite are linked

%% Import unique untargeted metabolites and parent compounds
clear;clc

mass_cutoff = 0.02;
rt_cutoff = 0.2;

target_compounds = readtable('targeted_compounds.csv','ReadRowNames',true);

load('novel_untargeted_metabolomics_table.mat');

drug_list = unique(unique_table.drug);

%Add new rows for analysis
col_list = {'link','metabolite_present','drug_present'};
unique_table{:,col_list} = ...
    nan(size(unique_table,1),length(col_list));

% Match drugs to metabolites and check for links
for i = 1:length(drug_list)
    drug = drug_list{i};
    drug_properties = target_compounds{drug,{'product_ion','RT'}};
    
    results_filename = ['GNPS_analysis/GNPS_results/',drug,'_1000.txt'];
    cluster_filename = ['GNPS_analysis/GNPS_results/',drug,'_cluster.txt'];
    drug_metabolite_index = find(contains(unique_table.drug,drug));
    GNPS_results = readtable(results_filename,'Delimiter','\t');
    cluster_info = readtable(cluster_filename,'Delimiter','\t');
    
    for k = 1:length(drug_metabolite_index)
        
        metabolite_properties_H = unique_table{drug_metabolite_index(k),{'m_plus_H','RT'}};
        metabolite_properties_Na = unique_table{drug_metabolite_index(k),{'m_plus_Na','RT'}};
        metabolite_properties_K = unique_table{drug_metabolite_index(k),{'m_plus_K','RT'}};
        
        [link_H,metabolite_present_H,drug_present] = ...
            check_metabolite_relatedness(drug_properties,metabolite_properties_H,...
            GNPS_results,cluster_info,mass_cutoff,rt_cutoff);
        
        [link_Na,metabolite_present_Na,~] = ...
            check_metabolite_relatedness(drug_properties,metabolite_properties_Na,...
            GNPS_results,cluster_info,mass_cutoff,rt_cutoff);
        
        [link_K,metabolite_present_K,~] = ...
            check_metabolite_relatedness(drug_properties,metabolite_properties_K,...
            GNPS_results,cluster_info,mass_cutoff,rt_cutoff);
        
        unique_table{drug_metabolite_index(k),'link'} = ...
            link_H | link_Na | link_K;
        unique_table{drug_metabolite_index(k),'metabolite_present'} = ...
            metabolite_present_H | metabolite_present_Na | metabolite_present_K;
        unique_table{drug_metabolite_index(k),'drug_present'} = drug_present;
    end
end



%% Export table with reformatted donor numbers and new column names
new_donors = arrayfun(@(x) num2str(x{1}'),unique_table.donors,'UniformOutput',false);

export_table = unique_table;
export_table.donors = new_donors;

export_table.Properties.VariableNames{1} = 'drug_name';
export_table.Properties.VariableNames{2} = 'metabolite_exact_mass_from_profinder';
export_table.Properties.VariableNames{3} = 'metabolite_RT';
export_table.Properties.VariableNames{4} = 'donors_producing_metabolite';
export_table.Properties.VariableNames{5} = 'calculated_metabolite_m_plus_H';
export_table.Properties.VariableNames{6} = 'calculated_metabolite_m_plus_Na';
export_table.Properties.VariableNames{7} = 'calculated_metabolite_m_plus_K';
export_table.Properties.VariableNames{8} = 'observed_drug_mz';
export_table.Properties.VariableNames{10} = 'donors_run_for_MolNet';
export_table.Properties.VariableNames{11} = 'drug_and_metabolite_clustered_in_MolNet';
export_table.Properties.VariableNames{12} = 'MolNet_metabolite_observed';
export_table.Properties.VariableNames{13} = 'MolNet_parent_drug_observed';

export_table.drug_and_metabolite_clustered_in_MolNet = ...
    num2cell(export_table.drug_and_metabolite_clustered_in_MolNet);
applicable_link = export_table.MolNet_metabolite_observed & ...
    export_table.MolNet_parent_drug_observed;
export_table.drug_and_metabolite_clustered_in_MolNet(~applicable_link) = {'NA'};

writetable(export_table,'unique_metabolites_with_networking.csv',...
    'Delimiter',',','WriteVariableNames',true);


