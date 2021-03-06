function [permutation_pval,parametric_pval_array] = ...
    compute_permutation_pvalue(combined,combos,grp,num_permutations,exp_tstat,true_combos)

%This function computes a permutation p value by shuffling sample labels. It
%uses a two-sided assumption. If the two classes being compared are
%identical (i.e. donors againsts donors), the script ensures that only
%non-redundant samples are counted.

if exist('true_combos')
    same_class = true;
else
    same_class = false;
end

tstat_array = zeros(num_permutations,1);
parametric_pval_array = zeros(num_permutations,1);

same_donor_combos = combos(grp,:);

old_samples1 = 1:max(combos(:,1));
old_samples2 = 1:max(combos(:,2));

for i = 1:num_permutations
    
    permuted_combos = zeros(size(combos));
    
    %Shuffle sample names
    new_samples1 = old_samples1(randperm(length(old_samples1)));
    new_samples2 = old_samples2(randperm(length(old_samples2)));
    
    %Get new combination identities
    permuted_combos(:,1) = arrayfun(@(x) new_samples1(x == old_samples1),combos(:,1));
    permuted_combos(:,2) = arrayfun(@(x) new_samples2(x == old_samples2),combos(:,2));
    
    %If samples are of the same class, remove redundant samples (reduce
    %combos to n choose k plus self comparisons)
    if same_class
        [permuted_combos,true_index] = intersect(permuted_combos,true_combos,'rows');
        new_combined = combined(true_index);
    else
        new_combined = combined;
    end
    
    %Figure out where the new same donor combos are
    [~,same_donor_indices] = intersect(permuted_combos,same_donor_combos,'rows');
    
    %Construct a new grp variable reflecting what is now a self-self combo
    permuted_grp = zeros(length(new_combined),1);
    permuted_grp(same_donor_indices) = 1;
    permuted_grp = logical(permuted_grp);
    
    %Compute parametric pval and get test stat
    [~,parametric_pval_array(i),~,stats] = ...
        ttest2(new_combined(permuted_grp), new_combined(~permuted_grp),'Vartype','unequal');
    
    tstat_array(i) = stats.tstat;
    
end

%Compute permutation p-val
permutation_pval = (1+sum(abs(tstat_array) > exp_tstat))/(1+num_permutations);

end

