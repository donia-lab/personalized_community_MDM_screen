function [permutation_pval,parametric_pval_array] = ...
    compute_permutation_pvalue(combined,combos,grp,num_permutations,exp_tstat)

%This function computes a permutation p value by shuffling sample labels. It
%uses a two-sided assumptions.

tstat_array = zeros(num_permutations,1);
parametric_pval_array = zeros(num_permutations,1);

same_donor_combos = combos(grp,:);

old_samples1 = 1:max(combos(:,1));
old_samples2 = 1:max(combos(:,2));

for i = 1:num_permutations
    
    new_samples1 = old_samples1(randperm(length(old_samples1)));
    new_samples2 = old_samples2(randperm(length(old_samples2)));
    
    permuted_combos(:,1) = arrayfun(@(x) new_samples1(x == old_samples1),combos(:,1));
    permuted_combos(:,2) = arrayfun(@(x) new_samples2(x == old_samples2),combos(:,2));
    
    [~,same_donor_indices] = intersect(permuted_combos,same_donor_combos,'rows');
    
    permuted_grp = zeros(length(combined),1);
    permuted_grp(same_donor_indices) = 1;
    permuted_grp = logical(permuted_grp);
    
    [~,parametric_pval_array(i),~,stats] = ...
        ttest2(combined(permuted_grp), combined(~permuted_grp),'Vartype','unequal');

    tstat_array(i) = stats.tstat;
    
end

permutation_pval = (1+sum(abs(tstat_array) > exp_tstat))/(1+num_permutations);

end

