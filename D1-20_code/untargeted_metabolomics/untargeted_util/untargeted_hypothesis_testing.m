function [BG_pval,HK_pval,DMSO_pval,BG_fc,HK_fc,DMSO_fc] = ...
    untargeted_hypothesis_testing(donor_manifest,BG_manifest,HK_manifest,compound_list)


%Enumerate all donor-compound combos for hypothesis testing
donors = unique(donor_manifest.donor);
compounds = unique(donor_manifest.compound);
[donor_combos,compound_combos] = enumerate_combos(donors,compounds);

%Make the big pval matrices
entity_num = size(compound_list,1);
entity_names = cellstr(strcat('entity_',string(1:entity_num)));
pval_varnames = ['donor','compound',entity_names];

BG_pval = array2table(NaN(length(donor_combos),length(pval_varnames)),...
    'VariableNames',pval_varnames);
BG_pval.donor = donor_combos;
BG_pval.compound = compound_combos;

HK_pval = BG_pval;
DMSO_pval = BG_pval;
BG_fc = BG_pval;
HK_fc = BG_pval;
DMSO_fc = BG_pval;

%Loop through and hypothesis test!

for i = 1:size(BG_pval,1)
    donor = BG_pval.donor(i);
    compound = BG_pval.compound{i};
    for j = 1:entity_num
        entity_name = ['entity_',num2str(j)];
        donor_replicates = extract_replicates(donor_manifest,donor,compound,j);
        BG_replicates = extract_replicates(BG_manifest,0,compound,j);
        HK_replicates = extract_replicates(HK_manifest,0,compound,j);
        DMSO_replicates = extract_replicates(donor_manifest,donor,'DMSO',j);
        
        BG_pval{i,entity_name} = wrap_ttest(donor_replicates,BG_replicates); 
        HK_pval{i,entity_name} = wrap_ttest(donor_replicates,HK_replicates); 
        DMSO_pval{i,entity_name} = wrap_ttest(donor_replicates,DMSO_replicates); 
        
        BG_fc{i,entity_name} = log2(mean(donor_replicates)/mean(BG_replicates));
        HK_fc{i,entity_name} = log2(mean(donor_replicates)/mean(HK_replicates));
        DMSO_fc{i,entity_name} = log2(mean(donor_replicates)/mean(DMSO_replicates));
    end
end

end

