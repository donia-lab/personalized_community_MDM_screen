function set_table = find_optimal_sample_set(metabolite_table)

%This script finds the minimal set of samples to be rerun to characterize
%all metabolites with molecular networking

drug_list = unique(metabolite_table.drug);
hitting_sets = {};
for i = 1:length(drug_list)
    drug = drug_list{i};
    drug_metabolites = ...
        metabolite_table(contains(metabolite_table.drug,drug),:); 
    donor_sets = drug_metabolites.donors;
    
    hitting_set = greedy_hitting_set(donor_sets');
    
    hitting_sets = [hitting_sets; hitting_set];
end

set_table = table(drug_list,hitting_sets);


end
