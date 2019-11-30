function [combined,grp,pval,exp_tstat] = ...
    compute_self_nonself_shared_ASVs(manifest1,manifest2,combos)

shared_entry = @(v1,v2) sum((v1 > 0)&(v2 > 0));

grp = zeros(size(combos,1),1);
combined = zeros(size(combos,1),1);

for i = 1:size(combos,1)
    i1 = combos(i,1);
    i2 = combos(i,2);
    v1 = manifest1.rel_asv{i1}.Variables;
    v1_donor = manifest1.donor(i1);
    v2 = manifest2.rel_asv{i2}.Variables;
    v2_donor = manifest2.donor(i2);
    
    grp(i) = v1_donor == v2_donor;
    combined(i) = shared_entry(v1,v2); 

end

grp = logical(grp);

[~,pval,~,stats] = ttest2(combined(grp), combined(~grp),'Vartype','unequal');

exp_tstat = stats.tstat;



end

